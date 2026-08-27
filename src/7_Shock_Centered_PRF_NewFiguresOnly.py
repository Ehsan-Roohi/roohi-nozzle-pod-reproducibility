# ============================================================
# 7_Shock_Centered_PRF_NewFiguresOnly.py
#
# Purpose
#   PRF-strengthening diagnostics ONLY.
#   This script intentionally DOES NOT regenerate the older Version-5 figures.
#
# Reads Tecplot ASCII files in the same folder as this script:
#   P=15.dat, P=20.dat, P=25.dat, P=30.dat, P=33.dat, ...
#
# Outputs are written only to:
#   ShockCenteredROM_PRF_V7/
#
# New outputs for the PRF revision:
#   1) Fig_PRF_V7_centerline_POD_full_spectrum_Rho.pdf/.png
#   2) Fig_PRF_V7_2D_POD_energy_Rho.pdf/.png
#   3) Fig_PRF_V7_2D_POD_modes_Rho.pdf/.png
#   4) Tab_PRF_V7_full_pod_spectrum_Rho.tex/.csv
#   5) Tab_PRF_V7_2D_POD_summary_Rho.tex/.csv
#   6) V7_run_summary.txt
#
# Main idea
#   - Reuse the original DSMC parser and shock-detection logic.
#   - Keep only the diagnostics that answer likely PRF reviewer concerns:
#       (i) full POD spectrum / number of modes for 95%, 99%, 99.5%,
#      (ii) 2D shock-window POD, not just centerline POD.
#
# Author note
#   Put this file in the folder containing P=15.dat, P=20.dat, ...
#   Then run:
#       python 7_Shock_Centered_PRF_NewFiguresOnly.py
# ============================================================

from pathlib import Path
import re
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d

# ============================================================
# Settings
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
OUTDIR = REPO_ROOT / "generated" / "ShockCenteredROM_PRF_V7"

# Pressures used in the PRF low-rank diagnostic.
# Change this list only if your paper uses a different selected set.
SELECTED_PRESSURES = ["15", "20", "25", "30", "33"]

SMOOTH_SIGMA = 0.7
DPI = 600
CENTERLINE_MODE = "max_y"  # "max_y" for symmetry centerline in half-domain; use "mid_y" if needed.

# If automatic shock detection catches a wrong peak, set e.g. (0.6e-4, 2.0e-4).
SHOCK_SEARCH_X_RANGE = None

# Reference length used by the Knudsen number stored in the Tecplot file.
# If DSMC Kn = lambda / Lref, use the same Lref here.
LREF_FOR_KNGLL = 5.7e-5

POD_FIELD = "Rho"

# Only the important PRF coordinates are retained here.
GRID_CONFIG = {
    "x_physical": {"min": 0.0, "max": 210.0, "n": 600, "xlabel": r"$x\;(\mu m)$"},
    "d_unscaled": {"min": -90.0, "max": 90.0, "n": 600, "xlabel": r"$x-x_s\;(\mu m)$"},
    "xi_jump": {"min": -8.0, "max": 8.0, "n": 600, "xlabel": r"$\xi_j=(x-x_s)/\delta_j$"},
    "xi_Lrho": {"min": -8.0, "max": 8.0, "n": 600, "xlabel": r"$\xi_{L\rho}=(x-x_s)/L_{\rho,s}$"},
}

COORD_LATEX = {
    "x_physical": r"\(x\)",
    "d_unscaled": r"\(x-x_s\)",
    "xi_jump": r"\((x-x_s)/\delta_j\)",
    "xi_Lrho": r"\((x-x_s)/L_{\rho,s}\)",
}

COORD_PLOT_LABEL = {
    "x_physical": r"$x$",
    "d_unscaled": r"$x-x_s$",
    "xi_jump": r"$\xi_j$",
    "xi_Lrho": r"$\xi_{L\rho}$",
}

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["mathtext.default"] = "regular"

# ============================================================
# File handling
# ============================================================

def is_clean_pressure_file(path: Path) -> bool:
    return re.fullmatch(r"P=\d+(\.\d+)?\.dat", path.name) is not None


def get_pressure_label(path: Path) -> str:
    m = re.fullmatch(r"P=([0-9]+(?:\.\d+)?)\.dat", Path(path).name)
    return m.group(1) if m else Path(path).stem


def get_pressure_value(path: Path) -> float:
    try:
        return float(get_pressure_label(path))
    except Exception:
        return 1.0e99

# ============================================================
# Tecplot parser
# ============================================================

def parse_variables(lines):
    variables_text = ""
    reading = False
    for line in lines:
        s = line.strip()
        if "VARIABLES" in s.upper():
            reading = True
            variables_text += " " + s
            continue
        if reading:
            if s.upper().startswith("ZONE"):
                break
            if re.match(r"^[\s\+\-]?\d", s):
                break
            variables_text += " " + s
    variables = re.findall(r'"([^"]+)"', variables_text)
    if not variables:
        variables = [
            "X", "Y", "Density", "QX", "QY", "T",
            "U", "V", "Txy", "Mach", "Pressure", "Kundsen"
        ]
    return variables


def clean_var_name(name):
    return name.lower().replace("_", "").replace(" ", "")


def find_var_index(variables, candidates, default_index=None, required=True):
    clean_vars = [clean_var_name(v) for v in variables]
    for cand in candidates:
        c = clean_var_name(cand)
        if c in clean_vars:
            return clean_vars.index(c)
    for cand in candidates:
        c = clean_var_name(cand)
        for i, v in enumerate(clean_vars):
            if c in v or v in c:
                return i
    if default_index is not None:
        return default_index
    if required:
        raise ValueError(f"Could not find variable {candidates}. Variables = {variables}")
    return None


def parse_zone_header(lines, start_idx):
    header_lines = []
    data_start = start_idx + 1
    for k in range(start_idx, min(start_idx + 20, len(lines))):
        s = lines[k].strip()
        if k > start_idx and re.match(r"^[\s\+\-]?\d", s):
            data_start = k
            break
        header_lines.append(s)
        if "DATAPACKING" in s.upper() or "DT=" in s.upper():
            data_start = k + 1
    header = " ".join(header_lines)
    i_match = re.search(r"\bI\s*=\s*(\d+)", header, flags=re.IGNORECASE)
    j_match = re.search(r"\bJ\s*=\s*(\d+)", header, flags=re.IGNORECASE)
    if i_match is None or j_match is None:
        raise ValueError(f"Could not parse I,J from zone header:\n{header}")
    I = int(i_match.group(1))
    J = int(j_match.group(1))
    name_match = re.search(r'T\s*=\s*"([^"]+)"', header, flags=re.IGNORECASE)
    zone_name = name_match.group(1) if name_match else "ZONE"
    return header, data_start, I, J, zone_name


def read_numeric_values(lines, start_idx, end_idx, n_values_expected):
    vals_all = []
    for line in lines[start_idx:end_idx]:
        s = line.strip()
        if not s:
            continue
        upper = s.upper()
        if upper.startswith("TITLE") or upper.startswith("VARIABLES"):
            continue
        if upper.startswith("ZONE"):
            break
        if "I=" in upper or "J=" in upper or "DATAPACKING" in upper or "DT=" in upper:
            continue
        if "VARSHARELIST" in upper:
            continue
        parts = s.replace(",", " ").split()
        try:
            vals = [float(p.replace("D", "E").replace("d", "e")) for p in parts]
        except ValueError:
            continue
        vals_all.extend(vals)
        if len(vals_all) >= n_values_expected:
            break
    vals_all = np.array(vals_all, dtype=float)
    if vals_all.size < n_values_expected:
        raise ValueError(f"Not enough numeric values. Expected {n_values_expected}, got {vals_all.size}")
    return vals_all[:n_values_expected]


def read_tecplot_zones(filename):
    filename = Path(filename)
    with open(filename, "r", errors="ignore") as f:
        lines = f.readlines()

    variables = parse_variables(lines)
    ix = find_var_index(variables, ["X"], default_index=0)
    iy = find_var_index(variables, ["Y"], default_index=1)
    irho = find_var_index(variables, ["Density", "rho", "Rho"], default_index=2)

    optional_candidates = {
        "U": ["U", "VelocityX", "Ux"],
        "V": ["V", "VelocityY", "Vy"],
        "T": ["T", "Temperature"],
        "Mach": ["Mach", "M"],
        "Pressure": ["Pressure", "P"],
        "Knudsen": ["Knudsen", "Kundsen", "KNUDSEN", "Kn", "K"],
    }
    optional = {k: find_var_index(variables, c, required=False) for k, c in optional_candidates.items()}

    zone_starts = [i for i, line in enumerate(lines) if line.strip().upper().startswith("ZONE")]
    if not zone_starts:
        raise ValueError(f"No ZONE found in {filename.name}")

    zones = []
    for z, start in enumerate(zone_starts):
        end = zone_starts[z + 1] if z + 1 < len(zone_starts) else len(lines)
        _, data_start, I, J, zone_name = parse_zone_header(lines, start)
        nvars = len(variables)
        npts = I * J
        raw = read_numeric_values(lines, data_start, end, npts * nvars)
        data = raw.reshape((npts, nvars))
        zone = {
            "name": zone_name,
            "I": I,
            "J": J,
            "X": data[:, ix].reshape(J, I),
            "Y": data[:, iy].reshape(J, I),
            "Rho": data[:, irho].reshape(J, I),
        }
        for key, idx in optional.items():
            if idx is not None:
                zone[key] = data[:, idx].reshape(J, I)
        zones.append(zone)
    return zones

# ============================================================
# Gradients and shock diagnostics
# ============================================================

def compute_curvilinear_gradient(X, Y, F, smooth_sigma=0.0):
    F2 = F.copy()
    if smooth_sigma and smooth_sigma > 0:
        F2 = gaussian_filter(F2, sigma=smooth_sigma)
    F_eta, F_xi = np.gradient(F2, edge_order=2)
    X_eta, X_xi = np.gradient(X, edge_order=2)
    Y_eta, Y_xi = np.gradient(Y, edge_order=2)
    J = X_xi * Y_eta - X_eta * Y_xi
    J_safe = J.copy()
    J_safe[np.abs(J_safe) < 1e-30] = np.nan
    dFdx = (F_xi * Y_eta - F_eta * Y_xi) / J_safe
    dFdy = (-F_xi * X_eta + F_eta * X_xi) / J_safe
    return dFdx, dFdy, np.sqrt(dFdx**2 + dFdy**2)


def add_diagnostics(zones):
    out = []
    eps = 1.0e-30
    for z in zones:
        zz = dict(z)
        dRdx, dRdy, grad_mag = compute_curvilinear_gradient(zz["X"], zz["Y"], zz["Rho"], SMOOTH_SIGMA)
        zz["dRdx"] = dRdx
        zz["dRdy"] = dRdy
        zz["grad_mag"] = grad_mag
        if "Knudsen" in zz:
            zz["KnGLL_rho"] = zz["Knudsen"] * LREF_FOR_KNGLL * grad_mag / (np.abs(zz["Rho"]) + eps)
        out.append(zz)
    return out


def choose_internal_zone(zones):
    # For your current files this selects the 100 x 60 internal-nozzle zone,
    # not the smaller external-plume zone.
    return max(zones, key=lambda z: z["I"] * z["J"])


def extract_centerline_profile(zone, field):
    X = zone["X"]
    Y = zone["Y"]
    F = zone[field]
    row_y = np.nanmedian(Y, axis=1)
    if CENTERLINE_MODE == "max_y":
        j = int(np.nanargmax(row_y))
    elif CENTERLINE_MODE == "mid_y":
        target = 0.5 * (np.nanmin(row_y) + np.nanmax(row_y))
        j = int(np.nanargmin(np.abs(row_y - target)))
    else:
        raise ValueError("CENTERLINE_MODE must be max_y or mid_y")
    x = X[j, :].copy()
    y = Y[j, :].copy()
    f = F[j, :].copy()
    order = np.argsort(x)
    return x[order], f[order], y[order]


def local_mean(x, q, x0, side, dx, nwin=10):
    if side == "left":
        mask = (x >= x0 - nwin * dx) & (x <= x0 - 3 * dx)
    else:
        mask = (x >= x0 + 3 * dx) & (x <= x0 + nwin * dx)
    vals = q[mask]
    vals = vals[np.isfinite(vals)]
    return float(np.nanmean(vals)) if vals.size else np.nan


def fwhm_width(x, g, idx_peak):
    g_abs = np.abs(g)
    peak = g_abs[idx_peak]
    if not np.isfinite(peak) or peak <= 0:
        return np.nan
    half = 0.5 * peak
    n = len(x)
    i_left = idx_peak
    while i_left > 0 and g_abs[i_left] >= half:
        i_left -= 1
    i_right = idx_peak
    while i_right < n - 1 and g_abs[i_right] >= half:
        i_right += 1
    if i_left == idx_peak or i_right == idx_peak:
        return np.nan
    return float(abs(x[i_right] - x[i_left]))


def gradient_moment_width(x, g, xs, window_factor=8.0, fallback_width=None):
    w = np.abs(g).copy()
    w[~np.isfinite(w)] = 0.0
    if fallback_width is None or not np.isfinite(fallback_width) or fallback_width <= 0:
        dx = np.nanmedian(np.diff(np.unique(x[np.isfinite(x)])))
        fallback_width = 10.0 * dx
    mask = np.isfinite(x) & (np.abs(x - xs) <= window_factor * fallback_width)
    if np.sum(mask) < 5 or np.nansum(w[mask]) <= 0:
        return np.nan
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    num = trapezoid(((x[mask] - xs) ** 2) * w[mask], x[mask])
    den = trapezoid(w[mask], x[mask])
    if den <= 0:
        return np.nan
    return float(np.sqrt(abs(num / den)))


def detect_shock_metrics(zone):
    x, rho, _ = extract_centerline_profile(zone, "Rho")
    _, drdx, _ = extract_centerline_profile(zone, "dRdx")
    finite = np.isfinite(x) & np.isfinite(rho) & np.isfinite(drdx)

    if SHOCK_SEARCH_X_RANGE is not None:
        xmin, xmax = SHOCK_SEARCH_X_RANGE
        finite &= (x >= xmin) & (x <= xmax)
    else:
        span = np.nanmax(x) - np.nanmin(x)
        finite &= (x >= np.nanmin(x) + 0.05 * span) & (x <= np.nanmin(x) + 0.95 * span)

    idx_candidates = np.where(finite)[0]
    if idx_candidates.size < 5:
        raise RuntimeError("Not enough valid points for shock detection.")

    idx_peak = idx_candidates[np.nanargmax(np.abs(drdx[finite]))]
    xs = float(x[idx_peak])
    gpeak = float(np.abs(drdx[idx_peak]))

    dx_arr = np.diff(np.unique(x[np.isfinite(x)]))
    dx_arr = dx_arr[dx_arr > 0]
    dx = float(np.nanmedian(dx_arr))

    rho_left = local_mean(x, rho, xs, "left", dx, nwin=12)
    rho_right = local_mean(x, rho, xs, "right", dx, nwin=12)
    delta_rho = abs(rho_left - rho_right)
    if not np.isfinite(delta_rho) or delta_rho <= 0:
        delta_rho = float(np.nanmax(rho) - np.nanmin(rho))

    delta_jump = delta_rho / gpeak if gpeak > 0 else np.nan
    delta_fwhm = fwhm_width(x, drdx, idx_peak)
    delta_moment = gradient_moment_width(x, drdx, xs, fallback_width=delta_fwhm)
    rho_s = float(rho[idx_peak])
    Lrho_s = abs(rho_s) / gpeak if gpeak > 0 else np.nan

    kn_s = np.nan
    lambda_s = np.nan
    kn_max = np.nan
    kngll_s = np.nan
    kngll_max = np.nan

    if "Knudsen" in zone:
        xk, kn, _ = extract_centerline_profile(zone, "Knudsen")
        kn_s = float(np.interp(xs, xk, kn))
        kn_max = float(np.nanmax(kn))
        lambda_s = kn_s * LREF_FOR_KNGLL

    if "KnGLL_rho" in zone:
        xgll, kngll, _ = extract_centerline_profile(zone, "KnGLL_rho")
        kngll_s = float(np.interp(xs, xgll, kngll))
        kngll_max = float(np.nanmax(kngll))

    # Fallbacks to keep the script stable.
    if not np.isfinite(delta_jump) or delta_jump <= 0:
        delta_jump = 5.0 * dx
    if not np.isfinite(delta_fwhm) or delta_fwhm <= 0:
        delta_fwhm = delta_jump
    if not np.isfinite(delta_moment) or delta_moment <= 0:
        delta_moment = delta_jump
    if not np.isfinite(Lrho_s) or Lrho_s <= 0:
        Lrho_s = delta_jump
    if not np.isfinite(lambda_s) or lambda_s <= 0:
        lambda_s = delta_jump

    return {
        "xs": xs,
        "gpeak": gpeak,
        "rho_left": rho_left,
        "rho_right": rho_right,
        "delta_rho": delta_rho,
        "delta_jump": float(delta_jump),
        "delta_fwhm": float(delta_fwhm),
        "delta_moment": float(delta_moment),
        "rho_s": rho_s,
        "Lrho_s": float(Lrho_s),
        "dx": dx,
        "kn_s": kn_s,
        "kn_max_centerline": kn_max,
        "lambda_s": float(lambda_s),
        "kngll_rho_s": kngll_s,
        "kngll_rho_max_centerline": kngll_max,
    }

# ============================================================
# Centerline profiles and 1D POD
# ============================================================

def normalize_jump(q, x, xs, dx):
    q_left = local_mean(x, q, xs, "left", dx, nwin=12)
    q_right = local_mean(x, q, xs, "right", dx, nwin=12)
    denom = q_left - q_right
    if not np.isfinite(denom) or abs(denom) < 1e-30:
        denom = np.nanmax(q) - np.nanmin(q)
    if not np.isfinite(denom) or abs(denom) < 1e-30:
        denom = 1.0
    return (q - q_right) / denom


def coordinate_values(x, metrics, coord_name):
    xs = metrics["xs"]
    if coord_name == "x_physical":
        return x * 1.0e6
    if coord_name == "d_unscaled":
        return (x - xs) * 1.0e6
    if coord_name == "xi_jump":
        return (x - xs) / metrics["delta_jump"]
    if coord_name == "xi_Lrho":
        return (x - xs) / metrics["Lrho_s"]
    raise ValueError(f"Unknown coordinate {coord_name}")


def collect_case(filename):
    p_label = get_pressure_label(filename)
    zones = add_diagnostics(read_tecplot_zones(filename))
    zone = choose_internal_zone(zones)
    metrics = detect_shock_metrics(zone)

    x, rho, _ = extract_centerline_profile(zone, POD_FIELD)
    rho_norm = normalize_jump(rho, x, metrics["xs"], metrics["dx"])

    return {
        "p_label": p_label,
        "p_value": float(p_label),
        "zone_name": zone["name"],
        "metrics": metrics,
        "centerline": {"x": x, "Rho": rho, "Rho_norm": rho_norm},
        "zone": zone,
    }


def build_centerline_matrix(cases, coord_name):
    cfg = GRID_CONFIG[coord_name]
    grid = np.linspace(cfg["min"], cfg["max"], cfg["n"])
    rows = []
    labels = []

    for case in cases:
        x = case["centerline"]["x"]
        qn = case["centerline"]["Rho_norm"]
        c = coordinate_values(x, case["metrics"], coord_name)

        mask = np.isfinite(c) & np.isfinite(qn)
        c = c[mask]
        qn = qn[mask]
        if c.size < 5:
            continue

        order = np.argsort(c)
        c = c[order]
        qn = qn[order]
        cu, idx = np.unique(c, return_index=True)
        qn = qn[idx]
        if cu.size < 5:
            continue

        f = interp1d(cu, qn, kind="linear", bounds_error=False, fill_value=np.nan)
        qi = f(grid)
        finite = np.isfinite(qi)
        if np.sum(finite) < 5:
            continue
        first = np.where(finite)[0][0]
        last = np.where(finite)[0][-1]
        qi[:first] = qi[first]
        qi[last + 1:] = qi[last]

        rows.append(qi)
        labels.append(case["p_label"])

    if not rows:
        raise RuntimeError(f"No valid centerline snapshots for coordinate {coord_name}")

    return grid, np.vstack(rows), labels


def pod_from_matrix(A):
    mean = np.nanmean(A, axis=0)
    A0 = A - mean
    A0 = np.nan_to_num(A0, nan=0.0)
    U, S, VT = np.linalg.svd(A0, full_matrices=False)
    energy = S**2
    energy_frac = energy / np.sum(energy) if np.sum(energy) > 0 else np.zeros_like(energy)
    cumulative = np.cumsum(energy_frac)
    return {"mean": mean, "S": S, "VT": VT, "energy_frac": energy_frac, "cumulative": cumulative}


def centerline_pod_metrics(cases, coord_name):
    grid, A, labels = build_centerline_matrix(cases, coord_name)
    pod = pod_from_matrix(A)
    pod["grid"] = grid
    pod["A"] = A
    pod["labels"] = labels
    return pod


def modes_needed_for_energy(cumulative, target):
    idx = np.where(cumulative >= target)[0]
    if idx.size == 0:
        return int(len(cumulative))
    return int(idx[0] + 1)

# ============================================================
# PRF output 1: full POD spectrum for centerline profiles
# ============================================================

def write_centerline_full_pod_spectrum_table(cases):
    rows = []
    pod_by_coord = {}

    for coord_name in GRID_CONFIG.keys():
        pod = centerline_pod_metrics(cases, coord_name)
        pod_by_coord[coord_name] = pod
        cum = pod["cumulative"]
        ef = pod["energy_frac"]
        rows.append({
            "coordinate": coord_name,
            "E1_percent": float(ef[0] * 100.0) if len(ef) else np.nan,
            "E1_plus_E2_percent": float(cum[1] * 100.0) if len(cum) > 1 else float(cum[0] * 100.0),
            "N95": modes_needed_for_energy(cum, 0.95),
            "N99": modes_needed_for_energy(cum, 0.99),
            "N995": modes_needed_for_energy(cum, 0.995),
        })

    csv_path = OUTDIR / "Tab_PRF_V7_full_pod_spectrum_Rho.csv"
    tex_path = OUTDIR / "Tab_PRF_V7_full_pod_spectrum_Rho.tex"

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["coordinate", "E1_percent", "E1_plus_E2_percent", "N95", "N99", "N995"]
        )
        writer.writeheader()
        writer.writerows(rows)

    tex_lines = []
    tex_lines.append(r"\begin{table}[t]")
    tex_lines.append(r"\centering")
    tex_lines.append(r"\caption{Full POD compactness of the centerline density profiles. \(N_{95}\), \(N_{99}\), and \(N_{99.5}\) denote the minimum number of modes required to capture 95\%, 99\%, and 99.5\% of the fluctuation energy, respectively.}")
    tex_lines.append(r"\label{tab:prf_v7_full_pod_spectrum}")
    tex_lines.append(r"\begin{tabular}{lccccc}")
    tex_lines.append(r"\toprule")
    tex_lines.append(r"Coordinate & \(E_1\) & \(E_{1+2}\) & \(N_{95}\) & \(N_{99}\) & \(N_{99.5}\) \\")
    tex_lines.append(r"\midrule")
    for r in rows:
        coord = COORD_LATEX[r["coordinate"]]
        tex_lines.append(
            f"{coord} & {r['E1_percent']:.2f}\\% & {r['E1_plus_E2_percent']:.2f}\\% & "
            f"{r['N95']} & {r['N99']} & {r['N995']} \\\\"
        )
    tex_lines.append(r"\bottomrule")
    tex_lines.append(r"\end{tabular}")
    tex_lines.append(r"\end{table}")
    tex_path.write_text("\n".join(tex_lines))

    # Figure: cumulative POD energy for all selected coordinates.
    fig, ax = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
    for coord_name, pod in pod_by_coord.items():
        modes = np.arange(1, len(pod["cumulative"]) + 1)
        ax.plot(modes, pod["cumulative"] * 100.0, "o-", lw=2, label=COORD_PLOT_LABEL[coord_name])
    ax.axhline(95, ls="--", lw=1.1)
    ax.axhline(99, ls="--", lw=1.1)
    ax.set_ylim(75, 101.5)
    ax.set_xlabel("Number of POD modes", fontweight="bold")
    ax.set_ylabel("Cumulative energy (%)", fontweight="bold")
    ax.set_title("Centerline density POD spectrum", fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=True, fontsize=9)
    fig.savefig(OUTDIR / "Fig_PRF_V7_centerline_POD_full_spectrum_Rho.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / "Fig_PRF_V7_centerline_POD_full_spectrum_Rho.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] wrote {csv_path.name}, {tex_path.name}")
    print("[OK] wrote Fig_PRF_V7_centerline_POD_full_spectrum_Rho.pdf/.png")
    return rows, pod_by_coord

# ============================================================
# PRF output 2: 2D shock-window POD
# ============================================================

def normalize_field_jump_2d(F, X, metrics):
    xs = metrics["xs"]
    dx = metrics["dx"]

    left_mask = (X < xs - 4.0 * dx) & (X > xs - 16.0 * dx)
    right_mask = (X > xs + 4.0 * dx) & (X < xs + 16.0 * dx)

    q_left = np.nanmean(F[left_mask]) if np.any(left_mask) else np.nan
    q_right = np.nanmean(F[right_mask]) if np.any(right_mask) else np.nan

    denom = q_left - q_right
    if not np.isfinite(denom) or abs(denom) < 1e-30:
        denom = np.nanpercentile(F, 95) - np.nanpercentile(F, 5)
    if not np.isfinite(denom) or abs(denom) < 1e-30:
        denom = 1.0

    return (F - q_right) / denom


def build_2d_shock_window_matrix(cases, field="Rho", xi_min=-6.0, xi_max=6.0, nxi=180, neta=80):
    """
    Build 2D snapshots in registered coordinates:
        xi_j = (x - xs) / delta_j
        eta  = normalized transverse coordinate.

    Interpolation is row-wise in xi, then column-wise in eta.
    This avoids scattered-grid dependencies and works for structured Tecplot zones.
    """
    xi_grid = np.linspace(xi_min, xi_max, nxi)
    eta_grid = np.linspace(0.0, 1.0, neta)

    snapshots = []
    labels = []

    for case in cases:
        zone = case["zone"]
        if field not in zone:
            continue

        X = zone["X"]
        Y = zone["Y"]
        F = zone[field]
        metrics = case["metrics"]

        xs = metrics["xs"]
        delta = metrics["delta_jump"]
        if not np.isfinite(delta) or delta <= 0:
            continue

        xi = (X - xs) / delta
        y_min = np.nanmin(Y)
        y_max = np.nanmax(Y)
        eta = (Y - y_min) / (y_max - y_min + 1.0e-300)
        Fn = normalize_field_jump_2d(F, X, metrics)

        row_stack = []
        eta_rows = []

        for j in range(Fn.shape[0]):
            xij = xi[j, :]
            qj = Fn[j, :]
            etaj = np.nanmean(eta[j, :])

            mask = np.isfinite(xij) & np.isfinite(qj)
            if np.sum(mask) < 5:
                continue

            xij = xij[mask]
            qj = qj[mask]
            order = np.argsort(xij)
            xij = xij[order]
            qj = qj[order]
            xiu, idx = np.unique(xij, return_index=True)
            qj = qj[idx]

            if xiu.size < 5:
                continue

            f_xi = interp1d(xiu, qj, kind="linear", bounds_error=False, fill_value=np.nan)
            qi = f_xi(xi_grid)
            finite = np.isfinite(qi)
            if np.sum(finite) < 5:
                continue

            first = np.where(finite)[0][0]
            last = np.where(finite)[0][-1]
            qi[:first] = qi[first]
            qi[last + 1:] = qi[last]

            row_stack.append(qi)
            eta_rows.append(etaj)

        if len(row_stack) < 5:
            continue

        row_stack = np.vstack(row_stack)
        eta_rows = np.array(eta_rows)

        Aeta = np.zeros((neta, nxi))
        order_eta = np.argsort(eta_rows)
        eta_sorted = eta_rows[order_eta]

        for k in range(nxi):
            qcol = row_stack[order_eta, k]
            eta_unique, idx = np.unique(eta_sorted, return_index=True)
            q_unique = qcol[idx]
            if eta_unique.size < 2:
                Aeta[:, k] = np.nanmean(q_unique)
            else:
                f_eta = interp1d(eta_unique, q_unique, kind="linear", bounds_error=False, fill_value="extrapolate")
                Aeta[:, k] = f_eta(eta_grid)

        snapshots.append(Aeta.reshape(-1))
        labels.append(case["p_label"])

    if not snapshots:
        raise RuntimeError("No valid 2D shock-window snapshots were built.")

    return xi_grid, eta_grid, np.vstack(snapshots), labels


def relative_reconstruction_errors(A, pod, ranks=(1, 2, 3)):
    mean = pod["mean"]
    VT = pod["VT"]
    A0 = np.nan_to_num(A - mean, nan=0.0)
    denom = np.linalg.norm(A0, axis=1)
    denom[denom < 1.0e-300] = 1.0

    out = {}
    for r in ranks:
        rr = min(r, VT.shape[0])
        Vr = VT[:rr, :]
        coeff = A0 @ Vr.T
        Arec = mean[None, :] + coeff @ Vr
        err = np.linalg.norm(A - Arec, axis=1) / denom
        out[r] = float(np.nanmean(err))
    return out


def write_2d_pod_results(cases, field="Rho"):
    xi_grid, eta_grid, A, labels = build_2d_shock_window_matrix(cases, field=field)
    pod = pod_from_matrix(A)

    ef = pod["energy_frac"]
    cum = pod["cumulative"]
    e1 = float(ef[0] * 100.0) if len(ef) else np.nan
    e12 = float(cum[1] * 100.0) if len(cum) > 1 else e1
    n95 = modes_needed_for_energy(cum, 0.95)
    n99 = modes_needed_for_energy(cum, 0.99)
    n995 = modes_needed_for_energy(cum, 0.995)
    rec = relative_reconstruction_errors(A, pod, ranks=(1, 2, 3))

    csv_path = OUTDIR / f"Tab_PRF_V7_2D_POD_summary_{field}.csv"
    tex_path = OUTDIR / f"Tab_PRF_V7_2D_POD_summary_{field}.tex"

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "field", "snapshots", "E1_percent", "E1_plus_E2_percent",
            "N95", "N99", "N995",
            "mean_rel_error_rank1", "mean_rel_error_rank2", "mean_rel_error_rank3"
        ])
        writer.writerow([field, len(labels), e1, e12, n95, n99, n995, rec[1], rec[2], rec[3]])

    tex_lines = []
    tex_lines.append(r"\begin{table}[t]")
    tex_lines.append(r"\centering")
    tex_lines.append(r"\caption{Two-dimensional shock-window POD of the density field in the registered coordinates \((\xi_j,\eta)\). Reconstruction errors are mean relative Frobenius errors over the DSMC snapshots.}")
    tex_lines.append(r"\label{tab:prf_v7_2d_pod_summary}")
    tex_lines.append(r"\begin{tabular}{ccccccc}")
    tex_lines.append(r"\toprule")
    tex_lines.append(r"\(E_1\) & \(E_{1+2}\) & \(N_{95}\) & \(N_{99}\) & \(N_{99.5}\) & rank-1 err. & rank-2 err. \\")
    tex_lines.append(r"\midrule")
    tex_lines.append(
        f"{e1:.2f}\\% & {e12:.2f}\\% & {n95} & {n99} & {n995} & {rec[1]:.3f} & {rec[2]:.3f} \\\\"
    )
    tex_lines.append(r"\bottomrule")
    tex_lines.append(r"\end{tabular}")
    tex_lines.append(r"\end{table}")
    tex_path.write_text("\n".join(tex_lines))

    # Figure 2D cumulative energy.
    fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
    modes = np.arange(1, len(cum) + 1)
    ax.plot(modes, cum * 100.0, "o-", lw=2)
    ax.axhline(95, ls="--", lw=1.1)
    ax.axhline(99, ls="--", lw=1.1)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Number of POD modes", fontweight="bold")
    ax.set_ylabel("Cumulative energy (%)", fontweight="bold")
    ax.set_title("2D shock-window density POD spectrum", fontweight="bold")
    ax.grid(True, alpha=0.3)
    fig.savefig(OUTDIR / f"Fig_PRF_V7_2D_POD_energy_{field}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / f"Fig_PRF_V7_2D_POD_energy_{field}.pdf", bbox_inches="tight")
    plt.close(fig)

    # Figure first three 2D modes in one compact article-ready figure.
    nxi = len(xi_grid)
    neta = len(eta_grid)
    nshow = min(3, pod["VT"].shape[0])
    fig, axes = plt.subplots(1, nshow, figsize=(4.2 * nshow, 3.5), constrained_layout=True)
    if nshow == 1:
        axes = [axes]
    for k, ax in enumerate(axes):
        mode = pod["VT"][k, :].reshape(neta, nxi)
        vmax = np.nanmax(np.abs(mode))
        levels = np.linspace(-vmax, vmax, 41) if vmax > 0 else 41
        im = ax.contourf(xi_grid, eta_grid, mode, levels=levels)
        ax.axvline(0.0, color="k", ls="--", lw=1.0)
        ax.set_xlabel(r"$\xi_j=(x-x_s)/\delta_j$", fontweight="bold")
        if k == 0:
            ax.set_ylabel(r"$\eta$", fontweight="bold")
        ax.set_title(rf"Mode {k+1}: {ef[k]*100.0:.2f}\%", fontweight="bold")
        fig.colorbar(im, ax=ax, shrink=0.86)
    fig.savefig(OUTDIR / f"Fig_PRF_V7_2D_POD_modes_{field}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / f"Fig_PRF_V7_2D_POD_modes_{field}.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] wrote {csv_path.name}, {tex_path.name}")
    print(f"[OK] wrote Fig_PRF_V7_2D_POD_energy_{field}.pdf/.png")
    print(f"[OK] wrote Fig_PRF_V7_2D_POD_modes_{field}.pdf/.png")

    return {
        "field": field,
        "labels": labels,
        "E1_percent": e1,
        "E1_plus_E2_percent": e12,
        "N95": n95,
        "N99": n99,
        "N995": n995,
        "rec": rec,
    }

# ============================================================
# Main
# ============================================================

def main():
    print("Working folder:")
    print(DATA_DIR)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    files = [f for f in DATA_DIR.glob("P=*.dat") if is_clean_pressure_file(f)]
    files = sorted(files, key=get_pressure_value)
    selected_files = [f for f in files if get_pressure_label(f) in SELECTED_PRESSURES]

    if not selected_files:
        raise FileNotFoundError(
            "No selected clean P=*.dat files found. Put this script in the folder with files such as P=15.dat, P=20.dat, ..."
        )

    print("\nSelected files:")
    for f in selected_files:
        print(" ", f.name)

    cases = []
    for f in selected_files:
        print(f"\nProcessing {f.name} ...")
        case = collect_case(f)
        cases.append(case)
        m = case["metrics"]
        print(
            f"  P={case['p_label']} kPa: xs={m['xs']*1e6:.3f} um, "
            f"delta_j={m['delta_jump']*1e6:.3f} um, "
            f"lambda_s={m['lambda_s']*1e6:.4f} um, "
            f"delta_j/lambda_s={m['delta_jump']/m['lambda_s']:.2f}"
        )

    spectrum_rows, _ = write_centerline_full_pod_spectrum_table(cases)
    summary_2d = write_2d_pod_results(cases, field=POD_FIELD)

    summary_path = OUTDIR / "V7_run_summary.txt"
    with open(summary_path, "w") as f:
        f.write("Version 7 PRF-only diagnostics\n")
        f.write("================================\n\n")
        f.write("Selected files:\n")
        for ff in selected_files:
            f.write(f"  {ff.name}\n")
        f.write("\nCenterline POD spectrum:\n")
        for r in spectrum_rows:
            f.write(
                f"  {r['coordinate']}: E1={r['E1_percent']:.2f}%, "
                f"E12={r['E1_plus_E2_percent']:.2f}%, "
                f"N95={r['N95']}, N99={r['N99']}, N99.5={r['N995']}\n"
            )
        f.write("\n2D shock-window POD:\n")
        f.write(
            f"  E1={summary_2d['E1_percent']:.2f}%, "
            f"E12={summary_2d['E1_plus_E2_percent']:.2f}%, "
            f"N95={summary_2d['N95']}, N99={summary_2d['N99']}, N99.5={summary_2d['N995']}\n"
        )
        f.write(
            f"  rank-1 error={summary_2d['rec'][1]:.4f}, "
            f"rank-2 error={summary_2d['rec'][2]:.4f}, "
            f"rank-3 error={summary_2d['rec'][3]:.4f}\n"
        )

    print(f"\n[OK] wrote {summary_path.name}")
    print("\nDone. Version-7 outputs only are in:")
    print(OUTDIR)


if __name__ == "__main__":
    main()

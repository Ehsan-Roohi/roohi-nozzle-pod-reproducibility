# ============================================================
# 5_Shock_Centered_ROM_MultiMetrics.py
#
# Purpose
#   Compare several reduced / shock-centered coordinates for DSMC
#   rarefied micro-nozzle centerline profiles.
#
# Reads Tecplot ASCII files:
#   P=15.dat, P=20.dat, P=25.dat, P=30.dat, P=33.dat, ...
#
# Outputs:
#   ShockCenteredROM_MultiMetrics/
#     shock_metrics_table.csv
#     coordinate_ranking.csv
#     coordinate_comparison_energy.png/pdf
#     coordinate_comparison_collapse_error.png/pdf
#     profiles_<coordinate>.png/pdf
#     pod_energy_<coordinate>_rho.png/pdf
#     pod_modes_<coordinate>_rho.png/pdf
#     kngll_centerline_<coordinate>.png/pdf  if Knudsen/Kundsen exists
#     prf_coordinate_compactness.png/pdf
#     prf_kinetic_scale_table.csv / .tex
#     prf_coordinate_compactness_table.csv / .tex
#
# Main coordinates tested
#   x_physical:      x
#   d_unscaled:      x - xs
#   xi_jump:         (x - xs) / delta_jump, delta_jump = Delta rho / max|drho/dx|
#   xi_fwhm:         (x - xs) / delta_fwhm, FWHM of |drho/dx| peak
#   xi_grad_moment:  (x - xs) / delta_moment, gradient-weighted rms width
#   xi_Lrho:         (x - xs) / Lrho_s, Lrho_s = |rho_s| / max|drho/dx|
#   xi_lambda:       (x - xs) / lambda_s, lambda_s = Kn_s * Lref
#   xi_gll_rho:      (x - xs) / Lrho_s, plus KnGLL diagnostic
#                    KnGLL_rho = Kn * Lref * |grad rho| / |rho|
#
# Notes
#   - The code recognizes both Knudsen and the common typo Kundsen.
#   - For neural-operator use, xs for held-out cases must be predicted from
#     training cases only. This script is diagnostic: it extracts xs from each
#     DSMC field to test whether reduced coordinates exist.
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
OUTDIR = REPO_ROOT / "generated" / "ShockCenteredROM_MultiMetrics"

SELECTED_PRESSURES = ["15", "20", "25", "30", "33"]

SMOOTH_SIGMA = 0.7
DPI = 600
CENTERLINE_MODE = "max_y"  # "max_y" or "mid_y"

# Put None for automatic. If automatic catches a wrong gradient peak,
# set this manually, e.g. (0.6e-4, 2.0e-4).
SHOCK_SEARCH_X_RANGE = None

# Reference length used by the Knudsen number stored in the file.
# If your DSMC Kn = lambda / Lref, use the same Lref here.
# If unsure, use the inlet height or characteristic nozzle height.
LREF_FOR_KNGLL = 5.7e-5

# Common grids for each coordinate.
GRID_CONFIG = {
    "x_physical":      {"min": 0.0,   "max": 210.0, "n": 600, "xlabel": r"$x\;(\mu m)$"},
    "d_unscaled":      {"min": -90.0, "max": 90.0,  "n": 600, "xlabel": r"$x-x_s\;(\mu m)$"},
    "xi_jump":         {"min": -8.0,  "max": 8.0,   "n": 600, "xlabel": r"$\xi_j=(x-x_s)/\delta_j$"},
    "xi_fwhm":         {"min": -8.0,  "max": 8.0,   "n": 600, "xlabel": r"$\xi_{fwhm}=(x-x_s)/\delta_{fwhm}$"},
    "xi_grad_moment":  {"min": -8.0,  "max": 8.0,   "n": 600, "xlabel": r"$\xi_m=(x-x_s)/\delta_m$"},
    "xi_Lrho":         {"min": -8.0,  "max": 8.0,   "n": 600, "xlabel": r"$\xi_{L\rho}=(x-x_s)/L_{\rho,s}$"},
    "xi_lambda":       {"min": -80.0, "max": 80.0,  "n": 600, "xlabel": r"$\xi_{\lambda}=(x-x_s)/\lambda_s$"},
}

FIELDS_TO_PLOT = ["Rho", "U", "Mach", "Pressure"]
FIELD_LABELS = {
    "Rho": r"$\rho$",
    "U": r"$U$",
    "Mach": r"$M$",
    "Pressure": r"$P$",
}

POD_FIELD = "Rho"

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["mathtext.default"] = "regular"

# ============================================================
# File handling
# ============================================================

def is_clean_pressure_file(path):
    return re.fullmatch(r"P=\d+(\.\d+)?\.dat", path.name) is not None


def get_pressure_label(path):
    m = re.fullmatch(r"P=([0-9]+(?:\.\d+)?)\.dat", Path(path).name)
    return m.group(1) if m else Path(path).stem


def get_pressure_value(path):
    try:
        return float(get_pressure_label(path))
    except Exception:
        return 1e99

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
# Gradient and diagnostics
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
    eps = 1e-30
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
        fallback_width = 10 * dx
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

    # fallbacks to prevent division by invalid scales
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
# Profile normalization and coordinates
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
        return x * 1e6
    if coord_name == "d_unscaled":
        return (x - xs) * 1e6
    if coord_name == "xi_jump":
        return (x - xs) / metrics["delta_jump"]
    if coord_name == "xi_fwhm":
        return (x - xs) / metrics["delta_fwhm"]
    if coord_name == "xi_grad_moment":
        return (x - xs) / metrics["delta_moment"]
    if coord_name == "xi_Lrho":
        return (x - xs) / metrics["Lrho_s"]
    if coord_name == "xi_lambda":
        return (x - xs) / metrics["lambda_s"]
    raise ValueError(f"Unknown coordinate {coord_name}")


def collect_case_profiles(filename):
    p_label = get_pressure_label(filename)
    zones = add_diagnostics(read_tecplot_zones(filename))
    zone = choose_internal_zone(zones)
    metrics = detect_shock_metrics(zone)
    profiles = {}
    available = [f for f in FIELDS_TO_PLOT if f in zone]
    if "Knudsen" in zone:
        available.append("Knudsen")
    if "KnGLL_rho" in zone:
        available.append("KnGLL_rho")
    if "abs_dRdx" not in available:
        zone["abs_dRdx"] = np.abs(zone["dRdx"])
        available.append("abs_dRdx")

    for field in available:
        x, q, _ = extract_centerline_profile(zone, field)
        if field in ["Knudsen", "KnGLL_rho", "abs_dRdx"]:
            qnorm = q / np.nanmax(np.abs(q)) if np.nanmax(np.abs(q)) > 0 else q
        else:
            qnorm = normalize_jump(q, x, metrics["xs"], metrics["dx"])
        profiles[field] = {"x": x, "q": q, "qnorm": qnorm}
    return {
        "p_label": p_label,
        "p_value": float(p_label),
        "zone_name": zone["name"],
        "metrics": metrics,
        "profiles": profiles,
    }

# ============================================================
# Interpolation, POD, and collapse metrics
# ============================================================

def pressure_text(case):
    return rf"$P_{{back}}={case['p_label']}\,\mathrm{{kPa}}$"


def build_matrix(cases, coord_name, field=POD_FIELD):
    cfg = GRID_CONFIG[coord_name]
    grid = np.linspace(cfg["min"], cfg["max"], cfg["n"])
    rows = []
    labels = []
    for case in cases:
        if field not in case["profiles"]:
            continue
        prof = case["profiles"][field]
        x = prof["x"]
        c = coordinate_values(x, case["metrics"], coord_name)
        qn = prof["qnorm"]
        mask = np.isfinite(c) & np.isfinite(qn)
        c = c[mask]
        qn = qn[mask]
        if c.size < 5:
            continue
        order = np.argsort(c)
        c = c[order]
        qn = qn[order]
        # Remove duplicate coordinates if any.
        cu, idx = np.unique(c, return_index=True)
        qn = qn[idx]
        if cu.size < 5:
            continue
        f = interp1d(cu, qn, kind="linear", bounds_error=False, fill_value=np.nan)
        qi = f(grid)
        # Fill edge NaNs with nearest finite values.
        finite = np.isfinite(qi)
        if np.sum(finite) < 5:
            continue
        qi[:np.argmax(finite)] = qi[finite][0]
        qi[np.where(finite)[0][-1] + 1:] = qi[finite][-1]
        rows.append(qi)
        labels.append(case["p_label"])
    A = np.vstack(rows)
    return grid, A, labels


def pod_metrics(cases, coord_name, field=POD_FIELD):
    grid, A, labels = build_matrix(cases, coord_name, field)
    mean = np.nanmean(A, axis=0)
    A0 = A - mean
    A0 = np.nan_to_num(A0, nan=0.0)
    U, S, VT = np.linalg.svd(A0, full_matrices=False)
    energy = S**2
    if np.sum(energy) <= 0:
        energy_frac = np.zeros_like(energy)
    else:
        energy_frac = energy / np.sum(energy)
    cumulative = np.cumsum(energy_frac)
    # Collapse error: mean RMS scatter around mean normalized by mean RMS magnitude.
    scatter = np.sqrt(np.mean((A - mean) ** 2, axis=0))
    collapse_error = float(np.sqrt(np.nanmean(scatter**2)))
    return {
        "grid": grid,
        "A": A,
        "labels": labels,
        "mean": mean,
        "S": S,
        "VT": VT,
        "energy_frac": energy_frac,
        "cumulative": cumulative,
        "collapse_error": collapse_error,
    }

def coord_label(coord_name):
    labels = {
        "x_physical": r"$x\;(\mu\mathrm{m})$",
        "d_unscaled": r"$x-x_s\;(\mu\mathrm{m})$",
        "xi_jump": r"$\xi_j=(x-x_s)/\delta_j$",
        "xi_fwhm": r"$\xi_{\mathrm{FWHM}}=(x-x_s)/\delta_{\mathrm{FWHM}}$",
        "xi_grad_moment": r"$\xi_m=(x-x_s)/\delta_m$",
        "xi_Lrho": r"$\xi_{L\rho}=(x-x_s)/L_{\rho,s}$",
        "xi_lambda": r"$\xi_{\lambda}=(x-x_s)/\lambda_s$",
    }
    return labels.get(coord_name, coord_name)


def coord_short_title(coord_name):
    labels = {
        "x_physical": r"physical coordinate $x$",
        "d_unscaled": r"translated coordinate $x-x_s$",
        "xi_jump": r"jump-scaled coordinate $\xi_j$",
        "xi_fwhm": r"FWHM-scaled coordinate $\xi_{\mathrm{FWHM}}$",
        "xi_grad_moment": r"gradient-moment coordinate $\xi_m$",
        "xi_Lrho": r"density-gradient-length coordinate $\xi_{L\rho}$",
        "xi_lambda": r"mean-free-path-scaled coordinate $\xi_{\lambda}$",
    }
    return labels.get(coord_name, coord_name)


def field_label(field):
    labels = {
        "Rho": r"density",
        "U": r"streamwise velocity",
        "Mach": r"Mach number",
        "Pressure": r"pressure",
    }
    return labels.get(field, field)

# ============================================================
# Plotting
# ============================================================

def plot_profiles_for_coordinate(cases, coord_name):
    cfg = GRID_CONFIG[coord_name]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    axes = axes.ravel()
    for ax, field in zip(axes, FIELDS_TO_PLOT):
        for case in cases:
            if field not in case["profiles"]:
                continue
            prof = case["profiles"][field]
            c = coordinate_values(prof["x"], case["metrics"], coord_name)
            ax.plot(c, prof["qnorm"], lw=2, label=pressure_text(case))
        if coord_name != "x_physical":
            ax.axvline(0.0, color="k", ls="--", lw=1.2)
        ax.set_xlim(cfg["min"], cfg["max"])
        ax.set_title(field, fontweight="bold")
        ax.set_xlabel(coord_label(coord_name), fontweight="bold")
        ax.set_ylabel(FIELD_LABELS.get(field, field) + " normalized", fontweight="bold")
        ax.grid(True, alpha=0.3)
    axes[0].legend(fontsize=9)
    fig.savefig(OUTDIR / f"profiles_{coord_name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / f"profiles_{coord_name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_pod_for_coordinate(pod, coord_name, field=POD_FIELD):
    cfg = GRID_CONFIG[coord_name]
    grid = pod["grid"]
    cumulative = pod["cumulative"]
    energy_frac = pod["energy_frac"]
    VT = pod["VT"]
    mean = pod["mean"]

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    modes = np.arange(1, len(cumulative) + 1)
    ax.plot(modes, cumulative * 100, "o-", lw=2)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Number of POD modes", fontweight="bold")
    ax.set_ylabel("Cumulative energy (%)", fontweight="bold")
    ax.set_title(
    r"POD energy of shock-centered density profiles",
    fontweight="bold",
    )
    ax.grid(True, alpha=0.3)
    fig.savefig(OUTDIR / f"pod_energy_{coord_name}_{field.lower()}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / f"pod_energy_{coord_name}_{field.lower()}.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.plot(grid, mean, lw=2.5, label="Mean profile")
    nshow = min(3, VT.shape[0])
    for k in range(nshow):
        ax.plot(grid, VT[k, :], lw=2, label=f"Mode {k+1}, energy={energy_frac[k]*100:.2f}%")
    if coord_name != "x_physical":
        ax.axvline(0.0, color="k", ls="--", lw=1.2)
    ax.set_xlim(cfg["min"], cfg["max"])
    ax.set_xlabel(coord_label(coord_name), fontweight="bold")
    ax.set_ylabel("POD mode amplitude", fontweight="bold")
    ax.set_title(
        rf"POD modes of shock-centered {field_label(field)} profiles",
        fontweight="bold",
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.savefig(OUTDIR / f"pod_modes_{coord_name}_{field.lower()}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / f"pod_modes_{coord_name}_{field.lower()}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_kngll_if_available(cases, coord_name):
    has = any("KnGLL_rho" in c["profiles"] for c in cases)
    if not has:
        return
    cfg = GRID_CONFIG[coord_name]
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    for case in cases:
        if "KnGLL_rho" not in case["profiles"]:
            continue
        prof = case["profiles"]["KnGLL_rho"]
        c = coordinate_values(prof["x"], case["metrics"], coord_name)
        ax.plot(c, prof["q"], lw=2, label=pressure_text(case))
    if coord_name != "x_physical":
        ax.axvline(0.0, color="k", ls="--", lw=1.2)
    ax.set_xlim(cfg["min"], cfg["max"])
    ax.set_xlabel(coord_label(coord_name), fontweight="bold")
    ax.set_ylabel(r"$Kn_{GLL,\rho}$", fontweight="bold")
    ax.set_title(
    rf"Gradient-length Knudsen number in {coord_short_title(coord_name)}",
    fontweight="bold",
    )
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.savefig(OUTDIR / f"kngll_centerline_{coord_name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / f"kngll_centerline_{coord_name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_coordinate_comparison(summary):
    coords = [coord_label(r["coordinate"]) for r in summary]
    e1 = [r["pod_mode1_energy_percent"] for r in summary]
    e2 = [r["pod_mode2_cumulative_percent"] for r in summary]
    cerr = [r["collapse_error"] for r in summary]

    fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=True)
    x = np.arange(len(coords))
    ax.bar(x - 0.18, e1, width=0.36, label="Mode 1 energy")
    ax.bar(x + 0.18, e2, width=0.36, label="Modes 1-2 cumulative")
    ax.set_xticks(x)
    ax.set_xticklabels(coords, rotation=35, ha="right")
    ax.set_ylabel("POD energy (%)", fontweight="bold")
    ax.set_title("Coordinate comparison by POD compactness", fontweight="bold")
    ax.set_ylim(0, 105)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.savefig(OUTDIR / "coordinate_comparison_energy.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / "coordinate_comparison_energy.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=True)
    ax.bar(x, cerr)
    ax.set_xticks(x)
    ax.set_xticklabels(coords, rotation=35, ha="right")
    ax.set_ylabel("Collapse error", fontweight="bold")
    ax.set_title("Coordinate comparison by profile scatter", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    fig.savefig(OUTDIR / "coordinate_comparison_collapse_error.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / "coordinate_comparison_collapse_error.pdf", bbox_inches="tight")
    plt.close(fig)

# ============================================================
# CSV outputs
# ============================================================

def write_metrics_table(cases):
    path = OUTDIR / "shock_metrics_table.csv"
    cols = [
        "P_back_kPa", "zone_name", "x_s_m", "x_s_um", "max_abs_drho_dx",
        "rho_left", "rho_right", "delta_rho", "delta_jump_m", "delta_jump_um",
        "delta_fwhm_m", "delta_fwhm_um", "delta_moment_m", "delta_moment_um",
        "rho_s", "Lrho_s_m", "Lrho_s_um", "Kn_s", "Kn_max_centerline",
        "lambda_s_m", "lambda_s_um", "delta_jump_over_lambda_s", "Lrho_s_over_lambda_s",
        "KnGLL_rho_s", "KnGLL_rho_max_centerline",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(cols)
        for c in cases:
            m = c["metrics"]
            writer.writerow([
                c["p_label"], c["zone_name"], m["xs"], m["xs"] * 1e6, m["gpeak"],
                m["rho_left"], m["rho_right"], m["delta_rho"], m["delta_jump"], m["delta_jump"] * 1e6,
                m["delta_fwhm"], m["delta_fwhm"] * 1e6, m["delta_moment"], m["delta_moment"] * 1e6,
                m["rho_s"], m["Lrho_s"], m["Lrho_s"] * 1e6, m["kn_s"], m["kn_max_centerline"],
                m["lambda_s"], m["lambda_s"] * 1e6,
                safe_div(m["delta_jump"], m["lambda_s"]),
                safe_div(m["Lrho_s"], m["lambda_s"]),
                m["kngll_rho_s"], m["kngll_rho_max_centerline"],
            ])
    print(f"[OK] wrote {path}")


def write_coordinate_ranking(summary):
    path = OUTDIR / "coordinate_ranking.csv"
    cols = ["coordinate", "pod_mode1_energy_percent", "pod_mode2_cumulative_percent", "collapse_error"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        for row in summary:
            writer.writerow({k: row[k] for k in cols})
    print(f"[OK] wrote {path}")


def safe_div(num, den):
    """Safe division for diagnostic ratios."""
    if den is None or (not np.isfinite(den)) or abs(den) < 1e-300:
        return np.nan
    return num / den


def write_prf_kinetic_scale_table(cases):
    """
    Write a compact PRF-oriented table connecting the shock-centered thickness
    scale to kinetic rarefaction measures.
    """
    csv_path = OUTDIR / "prf_kinetic_scale_table.csv"
    tex_path = OUTDIR / "prf_kinetic_scale_table.tex"

    rows = []
    for case in cases:
        m = case["metrics"]
        delta_over_lambda = safe_div(m["delta_jump"], m["lambda_s"])
        lrho_over_lambda = safe_div(m["Lrho_s"], m["lambda_s"])
        rows.append({
            "P_back_kPa": case["p_label"],
            "x_s_um": m["xs"] * 1e6,
            "delta_j_um": m["delta_jump"] * 1e6,
            "lambda_s_um": m["lambda_s"] * 1e6,
            "delta_j_over_lambda_s": delta_over_lambda,
            "Lrho_s_over_lambda_s": lrho_over_lambda,
            "Kn_s": m["kn_s"],
            "KnGLL_rho_s": m["kngll_rho_s"],
            "KnGLL_rho_max_centerline": m["kngll_rho_max_centerline"],
        })

    cols = [
        "P_back_kPa",
        "x_s_um",
        "delta_j_um",
        "lambda_s_um",
        "delta_j_over_lambda_s",
        "Lrho_s_over_lambda_s",
        "Kn_s",
        "KnGLL_rho_s",
        "KnGLL_rho_max_centerline",
    ]

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)

    def fmt(x, kind="float"):
        if x is None or not np.isfinite(float(x)):
            return "--"
        x = float(x)
        if kind == "int":
            return f"{x:.0f}"
        if abs(x) < 1e-2 or abs(x) >= 1e3:
            return f"{x:.3e}"
        return f"{x:.3f}"

    tex_lines = []
    tex_lines.append(r"\begin{table}[t]")
    tex_lines.append(r"\centering")
    tex_lines.append(r"\caption{Kinetic interpretation of the shock-centered thickness scale. Here \(\lambda_s=Kn_s L_{\mathrm{ref}}\) is the local mean free path evaluated at the compression-layer station, and \(Kn_{\mathrm{GLL},\rho}\) is the density-based gradient-length local Knudsen number.}")
    tex_lines.append(r"\label{tab:prf_kinetic_scale}")
    tex_lines.append(r"\begin{tabular}{ccccccc}")
    tex_lines.append(r"\toprule")
    tex_lines.append(r"\(P_{\mathrm{back}}\) (kPa) & \(x_s\) (\(\mu\)m) & \(\delta_j\) (\(\mu\)m) & \(\lambda_s\) (\(\mu\)m) & \(\delta_j/\lambda_s\) & \(Kn_s\) & \(Kn_{\mathrm{GLL},\rho}^{\max}\) \\")
    tex_lines.append(r"\midrule")
    for r in rows:
        tex_lines.append(
            f"{r['P_back_kPa']} & "
            f"{fmt(r['x_s_um'])} & "
            f"{fmt(r['delta_j_um'])} & "
            f"{fmt(r['lambda_s_um'])} & "
            f"{fmt(r['delta_j_over_lambda_s'])} & "
            f"{fmt(r['Kn_s'])} & "
            f"{fmt(r['KnGLL_rho_max_centerline'])} \\\\"
        )
    tex_lines.append(r"\bottomrule")
    tex_lines.append(r"\end{tabular}")
    tex_lines.append(r"\end{table}")

    tex_path.write_text("\n".join(tex_lines))

    print(f"[OK] wrote {csv_path}")
    print(f"[OK] wrote {tex_path}")


def write_prf_coordinate_compactness_table(summary):
    """
    Write the short PRF main-text table with only x, x-xs, xi_j, and xi_Lrho.
    """
    keep = ["x_physical", "d_unscaled", "xi_jump", "xi_Lrho"]
    label_map = {
        "x_physical": r"\(x\)",
        "d_unscaled": r"\(x-x_s\)",
        "xi_jump": r"\((x-x_s)/\delta_j\)",
        "xi_Lrho": r"\((x-x_s)/L_{\rho,s}\)",
    }

    rows = [r for r in summary if r["coordinate"] in keep]
    rows = sorted(rows, key=lambda r: keep.index(r["coordinate"]))

    csv_path = OUTDIR / "prf_coordinate_compactness_table.csv"
    tex_path = OUTDIR / "prf_coordinate_compactness_table.tex"

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["coordinate", "E1_percent", "E1_plus_E2_percent"])
        for r in rows:
            writer.writerow([
                label_map[r["coordinate"]].replace("\\", ""),
                r["pod_mode1_energy_percent"],
                r["pod_mode2_cumulative_percent"],
            ])

    tex_lines = []
    tex_lines.append(r"\begin{table}[t]")
    tex_lines.append(r"\centering")
    tex_lines.append(r"\caption{Compactness of centerline density profiles under progressively aligned coordinates. \(E_1\) is the energy captured by the leading POD mode, and \(E_{1+2}\) is the cumulative energy captured by the first two modes.}")
    tex_lines.append(r"\label{tab:prf_coordinate_pod_comparison}")
    tex_lines.append(r"\begin{tabular}{lcc}")
    tex_lines.append(r"\toprule")
    tex_lines.append(r"Coordinate & \(E_1\) & \(E_{1+2}\) \\")
    tex_lines.append(r"\midrule")
    for r in rows:
        e1 = r["pod_mode1_energy_percent"]
        e12 = r["pod_mode2_cumulative_percent"]
        coord = label_map[r["coordinate"]]
        e1txt = rf"\textbf{{{e1:.2f}\%}}" if r["coordinate"] == "xi_jump" else f"{e1:.2f}\\%"
        e12txt = rf"\textbf{{{e12:.2f}\%}}" if r["coordinate"] == "xi_Lrho" else f"{e12:.2f}\\%"
        tex_lines.append(f"{coord} & {e1txt} & {e12txt} \\\\")
    tex_lines.append(r"\bottomrule")
    tex_lines.append(r"\end{tabular}")
    tex_lines.append(r"\end{table}")

    tex_path.write_text("\n".join(tex_lines))

    print(f"[OK] wrote {csv_path}")
    print(f"[OK] wrote {tex_path}")


def plot_prf_coordinate_compactness(summary):
    """
    PRF-style compact figure with only the physically important coordinate
    hierarchy: x, x-xs, and (x-xs)/delta_j.
    """
    keep = ["x_physical", "d_unscaled", "xi_jump"]
    label_map = {
        "x_physical": r"$x$",
        "d_unscaled": r"$x-x_s$",
        "xi_jump": r"$(x-x_s)/\delta_j$",
    }

    rows = [r for r in summary if r["coordinate"] in keep]
    rows = sorted(rows, key=lambda r: keep.index(r["coordinate"]))

    labels = [label_map[r["coordinate"]] for r in rows]
    e1 = [r["pod_mode1_energy_percent"] for r in rows]
    e12 = [r["pod_mode2_cumulative_percent"] for r in rows]

    fig, ax = plt.subplots(figsize=(5.3, 3.5), constrained_layout=True)
    x = np.arange(len(labels))
    ax.bar(x - 0.18, e1, width=0.36, label=r"$E_1$")
    ax.bar(x + 0.18, e12, width=0.36, label=r"$E_{1+2}$")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("POD energy (%)", fontweight="bold")
    ax.set_ylim(75, 101)
    ax.set_title("Compactness of shock-centered profiles", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(frameon=True)

    fig.savefig(OUTDIR / "prf_coordinate_compactness.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / "prf_coordinate_compactness.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] wrote {OUTDIR / 'prf_coordinate_compactness.pdf'}")

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
        raise FileNotFoundError("No selected clean P=*.dat files found.")

    print("\nSelected files:")
    for f in selected_files:
        print(" ", f.name)

    cases = []
    for f in selected_files:
        print(f"\nProcessing {f.name} ...")
        case = collect_case_profiles(f)
        cases.append(case)
        m = case["metrics"]
        print(
            f"  P={case['p_label']} kPa: xs={m['xs']*1e6:.3f} um, "
            f"delta_jump={m['delta_jump']*1e6:.3f} um, "
            f"delta_fwhm={m['delta_fwhm']*1e6:.3f} um, "
            f"delta_moment={m['delta_moment']*1e6:.3f} um, "
            f"Kn_s={m['kn_s']:.4e}, KnGLL_s={m['kngll_rho_s']:.4e}"
        )

    write_metrics_table(cases)
    write_prf_kinetic_scale_table(cases)

    summary = []
    for coord_name in GRID_CONFIG.keys():
        print(f"\nAnalyzing coordinate: {coord_name}")
        plot_profiles_for_coordinate(cases, coord_name)
        plot_kngll_if_available(cases, coord_name)
        pod = pod_metrics(cases, coord_name, field=POD_FIELD)
        plot_pod_for_coordinate(pod, coord_name, field=POD_FIELD)
        e1 = float(pod["energy_frac"][0] * 100) if len(pod["energy_frac"]) > 0 else 0.0
        e2 = float(pod["cumulative"][1] * 100) if len(pod["cumulative"]) > 1 else e1
        cerr = float(pod["collapse_error"])
        summary.append({
            "coordinate": coord_name,
            "pod_mode1_energy_percent": e1,
            "pod_mode2_cumulative_percent": e2,
            "collapse_error": cerr,
        })
        print(f"  Mode 1 energy: {e1:.3f}%")
        print(f"  Modes 1-2 cumulative: {e2:.3f}%")
        print(f"  Collapse error: {cerr:.5e}")

    # Rank: high mode-1 energy and low collapse error are both useful.
    summary_sorted = sorted(summary, key=lambda r: (-r["pod_mode1_energy_percent"], r["collapse_error"]))
    write_coordinate_ranking(summary_sorted)
    plot_coordinate_comparison(summary_sorted)
    write_prf_coordinate_compactness_table(summary_sorted)
    plot_prf_coordinate_compactness(summary_sorted)

    print("\nCoordinate ranking:")
    for row in summary_sorted:
        print(
            f"  {row['coordinate']:<16s}  "
            f"E1={row['pod_mode1_energy_percent']:.3f}%  "
            f"E12={row['pod_mode2_cumulative_percent']:.3f}%  "
            f"collapse_error={row['collapse_error']:.5e}"
        )

    print("\nDone.")
    print(f"Outputs are in: {OUTDIR}")


if __name__ == "__main__":
    main()

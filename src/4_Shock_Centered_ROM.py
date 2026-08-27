# ============================================================
# 4_Shock_Centered_ROM.py
#
# Shock-centered reduced-order analysis for DSMC micro-nozzle
#
# Reads:
#   P=15.dat, P=20.dat, P=25.dat, P=30.dat, P=33.dat
#
# Main idea:
#   1) compute d(rho)/dx on curvilinear grid
#   2) detect shock/compression-layer station xs from centerline max |d rho/dx|
#   3) estimate effective shock thickness delta_s = Delta rho / max |d rho/dx|
#   4) transform x -> xi = (x - xs) / delta_s
#   5) show collapse of profiles in shock-centered coordinates
#   6) perform POD/SVD on aligned density profiles
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
OUTDIR = REPO_ROOT / "generated" / "ShockCenteredROM"

SELECTED_PRESSURES = ["15", "20", "25", "30", "33"]

SMOOTH_SIGMA = 0.7

# Use max_y if your half-domain has the symmetry line at maximum y.
# Use mid_y if you want a middle row instead.
CENTERLINE_MODE = "max_y"

# Restrict shock search range if needed.
# Based on your nozzle figures, this usually avoids inlet/throat/outlet artifacts.
# Set to None for automatic full-range search.
SHOCK_SEARCH_X_RANGE = None
# Example:
# SHOCK_SEARCH_X_RANGE = (0.5e-4, 2.4e-4)

# Common shock-centered grid for aligned profiles
XI_MIN = -8.0
XI_MAX = 8.0
N_XI = 500

DPI = 600


# ============================================================
# File handling
# ============================================================

def is_clean_pressure_file(path):
    return re.fullmatch(r"P=\d+(\.\d+)?\.dat", path.name) is not None


def get_pressure_label(path):
    m = re.fullmatch(r"P=([0-9]+(?:\.\d+)?)\.dat", Path(path).name)
    if m:
        return m.group(1)
    return Path(path).stem


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
            "U", "V", "Txy", "Mach", "Pressure", "Knudsen"
        ]

    return variables


def clean_var_name(name):
    return name.lower().replace("_", "").replace(" ", "")


def find_var_index(variables, candidates, default_index=None):
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

    raise ValueError(f"Could not find variable {candidates}. Variables = {variables}")


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

        if upper.startswith("TITLE"):
            continue
        if upper.startswith("VARIABLES"):
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
        raise ValueError(
            f"Not enough numeric values. Expected {n_values_expected}, got {vals_all.size}"
        )

    return vals_all[:n_values_expected]


def read_tecplot_zones(filename):
    filename = Path(filename)

    with open(filename, "r", errors="ignore") as f:
        lines = f.readlines()

    variables = parse_variables(lines)
    clean_vars = [clean_var_name(v) for v in variables]

    ix = find_var_index(variables, ["X"], default_index=0)
    iy = find_var_index(variables, ["Y"], default_index=1)
    irho = find_var_index(variables, ["Density", "rho", "Rho"], default_index=2)

    optional = {}
    for key, candidates in {
        "U": ["U", "VelocityX", "Ux"],
        "V": ["V", "VelocityY", "Vy"],
        "T": ["T", "Temperature"],
        "Mach": ["Mach", "M"],
        "Pressure": ["Pressure", "P"],
        "Knudsen": ["Knudsen", "Kundsen", "KNUDSEN", "Kn", "K"],
    }.items():
        try:
            optional[key] = find_var_index(variables, candidates)
        except Exception:
            optional[key] = None

    zone_starts = []
    for i, line in enumerate(lines):
        if line.strip().upper().startswith("ZONE"):
            zone_starts.append(i)

    if not zone_starts:
        raise ValueError(f"No ZONE found in {filename.name}")

    zones = []

    for z, start in enumerate(zone_starts):
        end = zone_starts[z + 1] if z + 1 < len(zone_starts) else len(lines)

        header, data_start, I, J, zone_name = parse_zone_header(lines, start)

        nvars = len(variables)
        npts = I * J
        expected = npts * nvars

        raw = read_numeric_values(lines, data_start, end, expected)
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
# Gradient calculation
# ============================================================

def compute_curvilinear_gradient(X, Y, F, smooth_sigma=0.0):
    F2 = F.copy()

    if smooth_sigma is not None and smooth_sigma > 0:
        F2 = gaussian_filter(F2, sigma=smooth_sigma)

    F_eta, F_xi = np.gradient(F2, edge_order=2)
    X_eta, X_xi = np.gradient(X, edge_order=2)
    Y_eta, Y_xi = np.gradient(Y, edge_order=2)

    J = X_xi * Y_eta - X_eta * Y_xi
    J_safe = J.copy()
    J_safe[np.abs(J_safe) < 1e-30] = np.nan

    dFdx = (F_xi * Y_eta - F_eta * Y_xi) / J_safe
    dFdy = (-F_xi * X_eta + F_eta * X_xi) / J_safe
    grad_mag = np.sqrt(dFdx**2 + dFdy**2)

    return dFdx, dFdy, grad_mag


def add_density_gradients(zones):
    out = []

    for z in zones:
        dRdx, dRdy, grad_mag = compute_curvilinear_gradient(
            z["X"], z["Y"], z["Rho"], smooth_sigma=SMOOTH_SIGMA
        )

        zz = dict(z)
        zz["dRdx"] = dRdx
        zz["dRdy"] = dRdy
        zz["grad_mag"] = grad_mag
        out.append(zz)

    return out


# ============================================================
# Shock-centered extraction
# ============================================================

def choose_internal_zone(zones):
    """
    Your nozzle file has an internal nozzle zone and a plume zone.
    Usually the internal nozzle zone has more cells, e.g. 100x60.
    """
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
        raise ValueError("CENTERLINE_MODE must be 'max_y' or 'mid_y'")

    x = X[j, :].copy()
    y = Y[j, :].copy()
    f = F[j, :].copy()

    order = np.argsort(x)

    return x[order], f[order], y[order]


def local_mean(x, q, x0, side, dx, nwin=8):
    if side == "left":
        mask = (x >= x0 - nwin * dx) & (x <= x0 - 3 * dx)
    else:
        mask = (x >= x0 + 3 * dx) & (x <= x0 + nwin * dx)

    vals = q[mask]
    vals = vals[np.isfinite(vals)]

    if vals.size < 2:
        return np.nan

    return float(np.nanmean(vals))


def detect_shock_metrics(zone):
    x, rho, _ = extract_centerline_profile(zone, "Rho")
    _, drdx, _ = extract_centerline_profile(zone, "dRdx")

    finite = np.isfinite(x) & np.isfinite(rho) & np.isfinite(drdx)

    if SHOCK_SEARCH_X_RANGE is not None:
        xmin, xmax = SHOCK_SEARCH_X_RANGE
        finite &= (x >= xmin) & (x <= xmax)
    else:
        # Avoid extreme end points.
        xlo = np.nanmin(x) + 0.05 * (np.nanmax(x) - np.nanmin(x))
        xhi = np.nanmin(x) + 0.95 * (np.nanmax(x) - np.nanmin(x))
        finite &= (x >= xlo) & (x <= xhi)

    if np.sum(finite) < 5:
        raise RuntimeError("Not enough valid points for shock detection.")

    idx_candidates = np.where(finite)[0]
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

    delta_s = delta_rho / gpeak if gpeak > 0 else np.nan

    if not np.isfinite(delta_s) or delta_s <= 0:
        delta_s = 5.0 * dx

    return {
        "xs": xs,
        "delta_s": float(delta_s),
        "gpeak": gpeak,
        "rho_left": rho_left,
        "rho_right": rho_right,
        "delta_rho": float(delta_rho),
        "dx": dx,
    }


def normalize_jump(q, x, xs, dx):
    q_left = local_mean(x, q, xs, "left", dx, nwin=12)
    q_right = local_mean(x, q, xs, "right", dx, nwin=12)

    denom = q_left - q_right

    if not np.isfinite(denom) or abs(denom) < 1e-30:
        denom = np.nanmax(q) - np.nanmin(q)

    if not np.isfinite(denom) or abs(denom) < 1e-30:
        denom = 1.0

    return (q - q_right) / denom


def collect_case_profiles(filename):
    p_label = get_pressure_label(filename)
    p_value = float(p_label)

    zones = read_tecplot_zones(filename)
    zones = add_density_gradients(zones)

    zone = choose_internal_zone(zones)
    metrics = detect_shock_metrics(zone)

    xs = metrics["xs"]
    delta_s = metrics["delta_s"]
    dx = metrics["dx"]

    fields = ["Rho", "U", "Mach", "Pressure"]
    available_fields = [f for f in fields if f in zone]

    profiles = {}

    for field in available_fields:
        x, q, y = extract_centerline_profile(zone, field)
        xi = (x - xs) / delta_s
        qn = normalize_jump(q, x, xs, dx)

        profiles[field] = {
            "x": x,
            "xi": xi,
            "q": q,
            "qnorm": qn,
        }

    # gradient magnitude / dRdx diagnostic
    xg, qg, _ = extract_centerline_profile(zone, "dRdx")
    profiles["abs_dRdx"] = {
        "x": xg,
        "xi": (xg - xs) / delta_s,
        "q": np.abs(qg),
        "qnorm": np.abs(qg) / np.nanmax(np.abs(qg)),
    }

    # optional KnGLL-like diagnostic if Knudsen is available
    kn_max = np.nan
    if "Knudsen" in zone:
        xkn, kn, _ = extract_centerline_profile(zone, "Knudsen")
        kn_max = float(np.nanmax(kn))

    return {
        "p_label": p_label,
        "p_value": p_value,
        "zone_name": zone["name"],
        "metrics": metrics,
        "profiles": profiles,
        "kn_centerline_max": kn_max,
    }


# ============================================================
# Plotting
# ============================================================

def pressure_text(case):
    return rf"$P_{{back}}={case['p_label']}\,\mathrm{{kPa}}$"


def plot_physical_x(cases):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    axes = axes.ravel()

    fields = ["Rho", "U", "Mach", "Pressure"]
    labels = {
        "Rho": r"$\rho$",
        "U": r"$U$",
        "Mach": r"$M$",
        "Pressure": r"$P$",
    }

    for ax, field in zip(axes, fields):
        for case in cases:
            if field not in case["profiles"]:
                continue
            prof = case["profiles"][field]
            ax.plot(prof["x"] * 1e6, prof["qnorm"], lw=2, label=pressure_text(case))

        ax.set_xlabel(r"$x\;(\mu m)$", fontweight="bold")
        ax.set_ylabel(labels[field] + " normalized", fontweight="bold")
        ax.set_title(field, fontweight="bold")
        ax.grid(True, alpha=0.3)

    axes[0].legend(fontsize=9)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "collapse_x_physical.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / "collapse_x_physical.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_shock_centered_xi(cases):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    axes = axes.ravel()

    fields = ["Rho", "U", "Mach", "Pressure"]
    labels = {
        "Rho": r"$\rho$",
        "U": r"$U$",
        "Mach": r"$M$",
        "Pressure": r"$P$",
    }

    for ax, field in zip(axes, fields):
        for case in cases:
            if field not in case["profiles"]:
                continue
            prof = case["profiles"][field]
            ax.plot(prof["xi"], prof["qnorm"], lw=2, label=pressure_text(case))

        ax.axvline(0.0, color="k", lw=1.2, ls="--")
        ax.set_xlim(XI_MIN, XI_MAX)
        ax.set_xlabel(r"$\xi=(x-x_s)/\delta_s$", fontweight="bold")
        ax.set_ylabel(labels[field] + " normalized", fontweight="bold")
        ax.set_title(field, fontweight="bold")
        ax.grid(True, alpha=0.3)

    axes[0].legend(fontsize=9)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "collapse_xi_shock_centered.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / "collapse_xi_shock_centered.pdf", bbox_inches="tight")
    plt.close(fig)


def build_aligned_matrix(cases, field="Rho"):
    xi_common = np.linspace(XI_MIN, XI_MAX, N_XI)
    rows = []
    labels = []

    for case in cases:
        if field not in case["profiles"]:
            continue

        prof = case["profiles"][field]
        xi = prof["xi"]
        qn = prof["qnorm"]

        mask = np.isfinite(xi) & np.isfinite(qn)
        xi = xi[mask]
        qn = qn[mask]

        order = np.argsort(xi)
        xi = xi[order]
        qn = qn[order]

        # Need enough support on the common xi interval.
        if xi.min() > XI_MIN or xi.max() < XI_MAX:
            print(f"[WARN] Case P={case['p_label']} does not fully cover xi range. Still interpolating with edge fill.")

        f = interp1d(
            xi,
            qn,
            kind="linear",
            bounds_error=False,
            fill_value=(qn[0], qn[-1]),
        )

        rows.append(f(xi_common))
        labels.append(case["p_label"])

    A = np.vstack(rows)
    return xi_common, A, labels


def perform_pod(cases, field="Rho"):
    xi, A, labels = build_aligned_matrix(cases, field=field)

    mean_profile = np.mean(A, axis=0)
    A_fluc = A - mean_profile

    U, S, VT = np.linalg.svd(A_fluc, full_matrices=False)

    energy = S**2
    energy_frac = energy / np.sum(energy)
    cumulative = np.cumsum(energy_frac)

    # Energy plot
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    modes = np.arange(1, len(energy_frac) + 1)

    ax.plot(modes, cumulative * 100, "o-", lw=2)
    ax.set_xlabel("Number of POD modes", fontweight="bold")
    ax.set_ylabel("Cumulative energy (%)", fontweight="bold")
    ax.set_title(f"POD energy of shock-centered {field} profiles", fontweight="bold")
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)

    fig.savefig(OUTDIR / f"pod_energy_{field.lower()}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / f"pod_energy_{field.lower()}.pdf", bbox_inches="tight")
    plt.close(fig)

    # Mode plot
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.plot(xi, mean_profile, lw=2.5, label="Mean profile")

    nshow = min(3, VT.shape[0])
    for k in range(nshow):
        ax.plot(xi, VT[k, :], lw=2, label=f"Mode {k+1}, energy={energy_frac[k]*100:.2f}%")

    ax.axvline(0.0, color="k", lw=1.2, ls="--")
    ax.set_xlim(XI_MIN, XI_MAX)
    ax.set_xlabel(r"$\xi=(x-x_s)/\delta_s$", fontweight="bold")
    ax.set_ylabel("POD mode / mean", fontweight="bold")
    ax.set_title(f"POD modes of shock-centered {field} profiles", fontweight="bold")
    ax.legend(loc="upper right", frameon=True, fontsize=11, handlelength=2.4, borderpad=0.8)
    ax.grid(True, alpha=0.3)

    fig.savefig(OUTDIR / f"pod_modes_{field.lower()}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTDIR / f"pod_modes_{field.lower()}.pdf", bbox_inches="tight")
    plt.close(fig)

    return {
        "xi": xi,
        "A": A,
        "labels": labels,
        "energy_frac": energy_frac,
        "cumulative": cumulative,
        "singular_values": S,
    }


def write_shock_table(cases):
    OUTDIR.mkdir(parents=True, exist_ok=True)

    path = OUTDIR / "shock_layer_table.csv"

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "P_back_kPa",
            "zone_name",
            "x_s_m",
            "x_s_um",
            "delta_s_m",
            "delta_s_um",
            "max_abs_drho_dx",
            "rho_left",
            "rho_right",
            "delta_rho",
            "dx_m",
            "Kn_centerline_max_if_available",
        ])

        for case in cases:
            m = case["metrics"]

            writer.writerow([
                case["p_label"],
                case["zone_name"],
                m["xs"],
                m["xs"] * 1e6,
                m["delta_s"],
                m["delta_s"] * 1e6,
                m["gpeak"],
                m["rho_left"],
                m["rho_right"],
                m["delta_rho"],
                m["dx"],
                case["kn_centerline_max"],
            ])

    print(f"[OK] wrote {path}")


# ============================================================
# Main
# ============================================================

def main():
    print("Working folder:")
    print(DATA_DIR)

    files = [f for f in DATA_DIR.glob("P=*.dat") if is_clean_pressure_file(f)]
    files = sorted(files, key=get_pressure_value)

    if not files:
        raise FileNotFoundError(
            f"No clean raw files found in {DATA_DIR}. Expected files like P=15.dat"
        )

    selected_files = []
    for f in files:
        p = get_pressure_label(f)
        if p in SELECTED_PRESSURES:
            selected_files.append(f)

    if not selected_files:
        raise RuntimeError("No selected pressure files found.")

    print("\nSelected files:")
    for f in selected_files:
        print(" ", f.name)

    OUTDIR.mkdir(parents=True, exist_ok=True)

    cases = []

    for f in selected_files:
        print(f"\nProcessing {f.name} ...")
        case = collect_case_profiles(f)
        cases.append(case)

        m = case["metrics"]
        print(
            f"  P={case['p_label']} kPa: "
            f"x_s={m['xs']*1e6:.3f} um, "
            f"delta_s={m['delta_s']*1e6:.3f} um, "
            f"max|drho/dx|={m['gpeak']:.4e}"
        )

    write_shock_table(cases)

    print("\nMaking physical-x comparison ...")
    plot_physical_x(cases)

    print("Making shock-centered xi collapse plot ...")
    plot_shock_centered_xi(cases)

    print("Performing POD/SVD on shock-centered density profiles ...")
    pod = perform_pod(cases, field="Rho")

    print("\nPOD cumulative energy:")
    for i, c in enumerate(pod["cumulative"][:5], start=1):
        print(f"  first {i} mode(s): {100*c:.3f}%")

    print("\nDone.")
    print(f"Outputs are in: {OUTDIR}")


if __name__ == "__main__":
    main()

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from data_manip.extraction.telemac_file import TelemacFile


# ============================================================
# B2 Tingsanchali - tau_c sensitivity analysis
#
# Compare:
#   tau_c = 0.001 Pa
#   tau_c = 0.010 Pa  (validated reference case)
#   tau_c = 0.100 Pa
#
# Fixed:
#   porosity = 0.40
#   k_d      = 5.0e-4 m/(Pa s)
#   n        = 0.018
#   Ks       = 55.55555556
#
# Outputs:
#   1) Numerical summary at 0, 5, 10, 20, 30 s
#   2) Maximum bed lowering versus time
#   3) Final longitudinal bed profiles
#   4) Reservoir water level versus time
#   5) Crest discharge versus time
#
# IMPORTANT:
#   This script does NOT modify or rerun TELEMAC.
#   It only reads the three completed SLF files.
# ============================================================


# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

ROOT = (
    r"F:\Telemac_Projects\Tuanzho_New"
    r"\16_B2_Tingsanchali_erosion_validation"
)

RUN_DIR = os.path.join(ROOT, "05_run")
FIG_DIR = os.path.join(ROOT, "08_figures")

os.makedirs(FIG_DIR, exist_ok=True)


GEO_FILE = os.path.join(
    RUN_DIR,
    "B2_Tingsanchali_geo.slf"
)


CASES = {
    "0.001": os.path.join(
        RUN_DIR,
        "B2_Tingsanchali_erosion_tauc_0p001_res.slf"
    ),
    "0.010": os.path.join(
        RUN_DIR,
        "B2_Tingsanchali_erosion_ref_res.slf"
    ),
    "0.100": os.path.join(
        RUN_DIR,
        "B2_Tingsanchali_erosion_tauc_0p100_res.slf"
    ),
}


TARGET_TIMES = [
    0.0,
    5.0,
    10.0,
    20.0,
    30.0,
]


# ------------------------------------------------------------
# PHYSICAL / ANALYSIS SETTINGS
# ------------------------------------------------------------

# Reservoir-level sampling region.
#
# Keep away from the dam itself so that the extracted value
# represents upstream reservoir level rather than the strongly
# disturbed overtopping zone.
#
# Original dam front crest edge: x = 0
# Upstream toe: x = -2.40 m
#
# Use:
#   -15 <= x <= -5 m
# over the full flume width.

RES_X_MIN = -15.0
RES_X_MAX = -5.0


# Crest discharge section:
#
# Vertical section at x = 0, spanning the complete 1 m width.

CREST_X = 0.0
Y_MIN = 0.0
Y_MAX = 1.0


# Longitudinal bed profile:
#
# Use the centreline y = 0.5 m.

PROFILE_Y = 0.5
PROFILE_X_MIN = -3.0
PROFILE_X_MAX = 5.0


# Wet/dry threshold used only for diagnostic purposes.

H_WET = 1.0e-4


# ------------------------------------------------------------
# BASIC CHECKS
# ------------------------------------------------------------

print()
print("=" * 78)
print("B2 TAU_C SENSITIVITY ANALYSIS")
print("=" * 78)

if not os.path.isfile(GEO_FILE):
    print("ERROR: Geometry file not found:")
    print(GEO_FILE)
    sys.exit(1)

for label, path in CASES.items():
    if not os.path.isfile(path):
        print("ERROR: Result file not found:")
        print(path)
        sys.exit(1)

print("All required SLF files found.")
print()


# ------------------------------------------------------------
# OPEN FILES
# ------------------------------------------------------------

geo = TelemacFile(GEO_FILE)

results = {}

for label, path in CASES.items():

    print(
        f"Opening tau_c = {label} Pa:"
    )

    print(
        " ",
        os.path.basename(path)
    )

    results[label] = TelemacFile(path)

print()


# ------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------

def nearest_record(res, target_time):
    """
    Return record index and actual TELEMAC time nearest to
    target_time.
    """

    times = np.asarray(
        res.times,
        dtype=float
    )

    idx = int(
        np.argmin(
            np.abs(
                times - target_time
            )
        )
    )

    return idx, float(times[idx])


def get_bottom(res, record):
    """
    Read BOTTOM at a given record.
    """

    return np.asarray(
        res.get_data_value(
            "BOTTOM",
            record
        ),
        dtype=float
    )


def get_depth(res, record):
    """
    Read WATER DEPTH at a given record.
    """

    return np.asarray(
        res.get_data_value(
            "WATER DEPTH",
            record
        ),
        dtype=float
    )


def get_free_surface(res, record):
    """
    Read FREE SURFACE at a given record.
    """

    return np.asarray(
        res.get_data_value(
            "FREE SURFACE",
            record
        ),
        dtype=float
    )


def reservoir_level(res, record):
    """
    Area-node mean free-surface elevation in a quiet
    upstream reservoir region.
    """

    x = np.asarray(
        res.tri.x,
        dtype=float
    )

    mask = (
        (x >= RES_X_MIN)
        & (x <= RES_X_MAX)
    )

    if not np.any(mask):
        raise RuntimeError(
            "No nodes found in reservoir sampling region."
        )

    eta = get_free_surface(
        res,
        record
    )

    return float(
        np.mean(
            eta[mask]
        )
    )


def max_bed_lowering(
    z_initial,
    z_current
):
    """
    Maximum positive bed lowering relative to initial bed.
    Numerical bed raising is ignored here and checked
    separately.
    """

    lowering = (
        z_initial
        - z_current
    )

    return float(
        np.max(lowering)
    )


def crest_discharge(res, record):
    """
    Integrate depth-averaged discharge through the vertical
    cross-section x = 0 from y = 0 to y = 1.

    TELEMAC get_timeseries_on_polyline returns:
        points,
        curvilinear_abscissa,
        values

    We extract:
        H
        U

    and integrate H*U along the section.

    Since the section is vertical and flow is primarily in
    the x direction, H*U is the normal discharge per unit
    width.
    """

    polyline = [
        [CREST_X, Y_MIN],
        [CREST_X, Y_MAX],
    ]

    poly_number = [
        100
    ]

    points_h, s_h, h_all = (
        res.get_timeseries_on_polyline(
            "WATER DEPTH",
            polyline,
            poly_number
        )
    )

    points_u, s_u, u_all = (
        res.get_timeseries_on_polyline(
            "VELOCITY U",
            polyline,
            poly_number
        )
    )

    s_h = np.asarray(
        s_h,
        dtype=float
    )

    s_u = np.asarray(
        s_u,
        dtype=float
    )

    h_all = np.asarray(
        h_all,
        dtype=float
    )

    u_all = np.asarray(
        u_all,
        dtype=float
    )

    if not np.allclose(
        s_h,
        s_u,
        rtol=0.0,
        atol=1.0e-12
    ):
        raise RuntimeError(
            "H and U polyline abscissas are inconsistent."
        )

    h = h_all[:, record]
    u = u_all[:, record]

    hu = h * u

    q = np.trapz(
        hu,
        s_h
    )

    return float(q)


def centreline_bottom_profile(
    res,
    record
):
    """
    Extract longitudinal BOTTOM profile along y = 0.5 m.
    """

    polyline = [
        [PROFILE_X_MIN, PROFILE_Y],
        [PROFILE_X_MAX, PROFILE_Y],
    ]

    poly_number = [
        320
    ]

    points, s, values = (
        res.get_timeseries_on_polyline(
            "BOTTOM",
            polyline,
            poly_number
        )
    )

    points = np.asarray(
        points,
        dtype=float
    )

    values = np.asarray(
        values,
        dtype=float
    )

    x = points[:, 0]

    z = values[:, record]

    return x, z


# ------------------------------------------------------------
# INITIAL BED
# ------------------------------------------------------------

z0_geo = np.asarray(
    geo.get_data_value(
        "BOTTOM",
        0
    ),
    dtype=float
)

print("=" * 78)
print("INITIAL-BED CONSISTENCY")
print("=" * 78)

for label, res in results.items():

    rec0, t0 = nearest_record(
        res,
        0.0
    )

    z0_case = get_bottom(
        res,
        rec0
    )

    diff = np.max(
        np.abs(
            z0_case
            - z0_geo
        )
    )

    print(
        f"tau_c = {label} Pa : "
        f"max |Z(t=0)-Zgeo| = "
        f"{diff:.6e} m"
    )

print()


# ------------------------------------------------------------
# NUMERICAL SUMMARY AT TARGET TIMES
# ------------------------------------------------------------

summary = {}

print("=" * 78)
print("TIME-SLICE SUMMARY")
print("=" * 78)

for label, res in results.items():

    summary[label] = {
        "time": [],
        "dzmax": [],
        "eta": [],
        "q": [],
        "raising": [],
        "below_zero": [],
    }

    print()
    print(
        f"tau_c = {label} Pa"
    )

    print(
        "-" * 78
    )

    for target_time in TARGET_TIMES:

        rec, actual_time = nearest_record(
            res,
            target_time
        )

        z = get_bottom(
            res,
            rec
        )

        dzmax = max_bed_lowering(
            z0_geo,
            z
        )

        eta = reservoir_level(
            res,
            rec
        )

        q = crest_discharge(
            res,
            rec
        )

        raising = float(
            np.max(
                z - z0_geo
            )
        )

        below_zero = float(
            np.min(z)
        )

        summary[label]["time"].append(
            actual_time
        )

        summary[label]["dzmax"].append(
            dzmax
        )

        summary[label]["eta"].append(
            eta
        )

        summary[label]["q"].append(
            q
        )

        summary[label]["raising"].append(
            raising
        )

        summary[label]["below_zero"].append(
            below_zero
        )

        print(
            f"t = {actual_time:6.2f} s | "
            f"max lowering = {dzmax:10.6f} m | "
            f"eta_res = {eta:10.6f} m | "
            f"Qcrest = {q:12.8f} m3/s"
        )

    max_raise = max(
        summary[label]["raising"]
    )

    min_bottom = min(
        summary[label]["below_zero"]
    )

    print()
    print(
        f"Maximum bed raising diagnostic = "
        f"{max_raise:.6e} m"
    )

    print(
        f"Minimum BOTTOM over sampled times = "
        f"{min_bottom:.6f} m"
    )


# ------------------------------------------------------------
# FULL-TIME SERIES
# ------------------------------------------------------------

full_series = {}

print()
print("=" * 78)
print("BUILDING FULL TIME SERIES")
print("=" * 78)

for label, res in results.items():

    times = np.asarray(
        res.times,
        dtype=float
    )

    dz_series = []
    eta_series = []
    q_series = []

    for rec in range(len(times)):

        z = get_bottom(
            res,
            rec
        )

        dz_series.append(
            max_bed_lowering(
                z0_geo,
                z
            )
        )

        eta_series.append(
            reservoir_level(
                res,
                rec
            )
        )

        q_series.append(
            crest_discharge(
                res,
                rec
            )
        )

    full_series[label] = {
        "time": times,
        "dzmax": np.asarray(
            dz_series,
            dtype=float
        ),
        "eta": np.asarray(
            eta_series,
            dtype=float
        ),
        "q": np.asarray(
            q_series,
            dtype=float
        ),
    }

    print(
        f"tau_c = {label} Pa : "
        f"{len(times)} records processed."
    )


# ------------------------------------------------------------
# FINAL DIFFERENCES RELATIVE TO REFERENCE
# ------------------------------------------------------------

REF_LABEL = "0.010"

ref_res = results[REF_LABEL]

ref_rec, ref_time = nearest_record(
    ref_res,
    30.0
)

ref_z30 = get_bottom(
    ref_res,
    ref_rec
)

ref_eta30 = summary[REF_LABEL]["eta"][-1]
ref_q30 = summary[REF_LABEL]["q"][-1]
ref_dz30 = summary[REF_LABEL]["dzmax"][-1]


print()
print("=" * 78)
print("30-S DIFFERENCES RELATIVE TO tau_c = 0.010 Pa")
print("=" * 78)

for label, res in results.items():

    rec, actual_time = nearest_record(
        res,
        30.0
    )

    z30 = get_bottom(
        res,
        rec
    )

    dz30 = summary[label]["dzmax"][-1]
    eta30 = summary[label]["eta"][-1]
    q30 = summary[label]["q"][-1]

    bed_abs_diff = np.abs(
        z30 - ref_z30
    )

    max_bed_diff = float(
        np.max(
            bed_abs_diff
        )
    )

    mean_bed_diff = float(
        np.mean(
            bed_abs_diff
        )
    )

    dz_diff = (
        dz30 - ref_dz30
    )

    eta_diff = (
        eta30 - ref_eta30
    )

    q_diff = (
        q30 - ref_q30
    )

    if abs(ref_dz30) > 1.0e-12:
        dz_pct = (
            100.0
            * dz_diff
            / ref_dz30
        )
    else:
        dz_pct = np.nan

    if abs(ref_q30) > 1.0e-12:
        q_pct = (
            100.0
            * q_diff
            / ref_q30
        )
    else:
        q_pct = np.nan

    print()
    print(
        f"tau_c = {label} Pa"
    )

    print(
        f"  max |Z-Zref|      = "
        f"{max_bed_diff:.6e} m"
    )

    print(
        f"  mean |Z-Zref|     = "
        f"{mean_bed_diff:.6e} m"
    )

    print(
        f"  delta max erosion = "
        f"{dz_diff:+.6e} m "
        f"({dz_pct:+.3f} %)"
    )

    print(
        f"  delta eta_res     = "
        f"{eta_diff:+.6e} m"
    )

    print(
        f"  delta Qcrest      = "
        f"{q_diff:+.6e} m3/s "
        f"({q_pct:+.3f} %)"
    )


# ------------------------------------------------------------
# FIGURE 1
# MAXIMUM BED LOWERING VS TIME
# ------------------------------------------------------------

plt.figure(
    figsize=(8.0, 5.5)
)

for label in [
    "0.001",
    "0.010",
    "0.100",
]:

    data = full_series[label]

    plt.plot(
        data["time"],
        data["dzmax"],
        linewidth=2.0,
        label=rf"$\tau_c={label}$ Pa"
    )

plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    "Maximum bed lowering (m)"
)

plt.title(
    "B-2 sensitivity of maximum erosion to critical shear stress"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()

fig1 = os.path.join(
    FIG_DIR,
    "B2_tauc_max_erosion_vs_time.png"
)

plt.savefig(
    fig1,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FIGURE 2
# FINAL BED PROFILES AT 30 s
# ------------------------------------------------------------

plt.figure(
    figsize=(9.0, 5.5)
)

# Initial profile from reference case at t = 0

rec0, _ = nearest_record(
    ref_res,
    0.0
)

x0, z_initial_profile = (
    centreline_bottom_profile(
        ref_res,
        rec0
    )
)

plt.plot(
    x0,
    z_initial_profile,
    linewidth=2.0,
    linestyle="--",
    label="Initial bed"
)


for label in [
    "0.001",
    "0.010",
    "0.100",
]:

    res = results[label]

    rec30, _ = nearest_record(
        res,
        30.0
    )

    x, z = (
        centreline_bottom_profile(
            res,
            rec30
        )
    )

    plt.plot(
        x,
        z,
        linewidth=2.0,
        label=rf"$\tau_c={label}$ Pa"
    )


plt.xlabel(
    "x (m)"
)

plt.ylabel(
    "Bed elevation (m)"
)

plt.title(
    "B-2 longitudinal bed profiles at 30 s"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()

fig2 = os.path.join(
    FIG_DIR,
    "B2_tauc_final_bed_profiles.png"
)

plt.savefig(
    fig2,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FIGURE 3
# RESERVOIR LEVEL
# ------------------------------------------------------------

plt.figure(
    figsize=(8.0, 5.5)
)

for label in [
    "0.001",
    "0.010",
    "0.100",
]:

    data = full_series[label]

    plt.plot(
        data["time"],
        data["eta"],
        linewidth=2.0,
        label=rf"$\tau_c={label}$ Pa"
    )


plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    "Reservoir water level (m)"
)

plt.title(
    "B-2 reservoir-level sensitivity to critical shear stress"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()

fig3 = os.path.join(
    FIG_DIR,
    "B2_tauc_reservoir_level.png"
)

plt.savefig(
    fig3,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FIGURE 4
# CREST DISCHARGE
# ------------------------------------------------------------

plt.figure(
    figsize=(8.0, 5.5)
)

for label in [
    "0.001",
    "0.010",
    "0.100",
]:

    data = full_series[label]

    plt.plot(
        data["time"],
        data["q"],
        linewidth=2.0,
        label=rf"$\tau_c={label}$ Pa"
    )


plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    r"Crest discharge $Q$ ($m^3/s$)"
)

plt.title(
    "B-2 crest-discharge sensitivity to critical shear stress"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()

fig4 = os.path.join(
    FIG_DIR,
    "B2_tauc_crest_discharge.png"
)

plt.savefig(
    fig4,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FINAL REPORT
# ------------------------------------------------------------

print()
print("=" * 78)
print("FIGURES SAVED")
print("=" * 78)

print(
    os.path.basename(fig1)
)

print(
    os.path.basename(fig2)
)

print(
    os.path.basename(fig3)
)

print(
    os.path.basename(fig4)
)

print()

print(
    "IMPORTANT:"
)

print(
    "This script quantifies sensitivity only."
)

print(
    "Do not declare tau_c insensitive until the numerical "
    "differences and plots have been reviewed."
)

print()

print(
    "STATUS: B2 TAU_C SENSITIVITY ANALYSIS COMPLETED"
)

print("=" * 78)


# ------------------------------------------------------------
# CLOSE FILES
# ------------------------------------------------------------

for res in results.values():
    res.close()

geo.close()
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from data_manip.extraction.telemac_file import TelemacFile


# ============================================================
# B2 Tingsanchali - unified k_d analysis
#
# Compare five k_d cases:
#
#   KD050 = 0.5e-4 m/(Pa s)
#   KD100 = 1.0e-4 m/(Pa s)
#   KD250 = 2.5e-4 m/(Pa s)
#   KD500 = 5.0e-4 m/(Pa s)  [existing reference]
#   KD800 = 8.0e-4 m/(Pa s)
#
# Fixed:
#   tau_c    = 0.01 Pa
#   porosity = 0.40
#   n        = 0.018
#   Ks       = 55.55555556
#
# Outputs:
#   1. Numerical summary at 0, 5, 10, 20, 30 s
#   2. Numerical physical-consistency checks
#   3. Maximum erosion versus time
#   4. Longitudinal bed profiles at 30 s
#   5. Reservoir water level versus time
#   6. Crest discharge versus time
#
# IMPORTANT:
#   This script only analyzes completed TELEMAC results.
#   It does NOT modify or rerun any case.
# ============================================================


# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

ROOT = (
    r"F:\Telemac_Projects\Tuanzho_New"
    r"\16_B2_Tingsanchali_erosion_validation"
)

RUN_DIR = os.path.join(
    ROOT,
    "05_run"
)

FIG_DIR = os.path.join(
    ROOT,
    "08_figures"
)

os.makedirs(
    FIG_DIR,
    exist_ok=True
)


GEO_FILE = os.path.join(
    RUN_DIR,
    "B2_Tingsanchali_geo.slf"
)


# ------------------------------------------------------------
# CASE DEFINITIONS
# ------------------------------------------------------------

CASES = {
    "KD050": {
        "kd": 0.5e-4,
        "label": r"$k_d=0.5\times10^{-4}$",
        "file": os.path.join(
            RUN_DIR,
            "B2_Tingsanchali_erosion_kd_0p5e4_res.slf"
        ),
    },

    "KD100": {
        "kd": 1.0e-4,
        "label": r"$k_d=1.0\times10^{-4}$",
        "file": os.path.join(
            RUN_DIR,
            "B2_Tingsanchali_erosion_kd_1p0e4_res.slf"
        ),
    },

    "KD250": {
        "kd": 2.5e-4,
        "label": r"$k_d=2.5\times10^{-4}$",
        "file": os.path.join(
            RUN_DIR,
            "B2_Tingsanchali_erosion_kd_2p5e4_res.slf"
        ),
    },

    "KD500": {
        "kd": 5.0e-4,
        "label": r"$k_d=5.0\times10^{-4}$",
        "file": os.path.join(
            RUN_DIR,
            "B2_Tingsanchali_erosion_ref_res.slf"
        ),
    },

    "KD800": {
        "kd": 8.0e-4,
        "label": r"$k_d=8.0\times10^{-4}$",
        "file": os.path.join(
            RUN_DIR,
            "B2_Tingsanchali_erosion_kd_8p0e4_res.slf"
        ),
    },
}


CASE_ORDER = [
    "KD050",
    "KD100",
    "KD250",
    "KD500",
    "KD800",
]


TARGET_TIMES = [
    0.0,
    5.0,
    10.0,
    20.0,
    30.0,
]


# ------------------------------------------------------------
# ANALYSIS REGIONS
# ------------------------------------------------------------

# Quiet upstream reservoir region
RES_X_MIN = -15.0
RES_X_MAX = -5.0

# Crest discharge section
CREST_X = 0.0
Y_MIN = 0.0
Y_MAX = 1.0

# Longitudinal centreline profile
PROFILE_Y = 0.5
PROFILE_X_MIN = -3.0
PROFILE_X_MAX = 5.0


# ------------------------------------------------------------
# FILE CHECK
# ------------------------------------------------------------

print()
print("=" * 82)
print("B2 UNIFIED K_D ANALYSIS")
print("=" * 82)


if not os.path.isfile(GEO_FILE):

    print("ERROR: Geometry file not found:")
    print(GEO_FILE)

    sys.exit(1)


for key in CASE_ORDER:

    path = CASES[key]["file"]

    if not os.path.isfile(path):

        print("ERROR: Result file not found:")
        print(path)

        sys.exit(1)


print("Geometry + all five k_d result files found.")
print()


# ------------------------------------------------------------
# OPEN FILES
# ------------------------------------------------------------

geo = TelemacFile(
    GEO_FILE
)

results = {}


for key in CASE_ORDER:

    case = CASES[key]

    print(
        f"Opening {key}: "
        f"k_d = {case['kd']:.1e}"
    )

    results[key] = TelemacFile(
        case["file"]
    )


print()


# ------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------

def nearest_record(
    res,
    target_time
):

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


def get_bottom(
    res,
    record
):

    return np.asarray(
        res.get_data_value(
            "BOTTOM",
            record
        ),
        dtype=float
    )


def get_free_surface(
    res,
    record
):

    return np.asarray(
        res.get_data_value(
            "FREE SURFACE",
            record
        ),
        dtype=float
    )


def reservoir_level(
    res,
    record
):

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

    lowering = (
        z_initial
        - z_current
    )

    return float(
        np.max(
            lowering
        )
    )


def crest_discharge(
    res,
    record
):

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
# INITIAL GEOMETRY
# ------------------------------------------------------------

z0_geo = np.asarray(
    geo.get_data_value(
        "BOTTOM",
        0
    ),
    dtype=float
)


print("=" * 82)
print("INITIAL-BED CONSISTENCY")
print("=" * 82)


initial_ok = True


for key in CASE_ORDER:

    res = results[key]

    rec0, t0 = nearest_record(
        res,
        0.0
    )

    z0_case = get_bottom(
        res,
        rec0
    )


    diff = float(
        np.max(
            np.abs(
                z0_case - z0_geo
            )
        )
    )


    print(
        f"{key:6s} | "
        f"max |Z(t=0)-Zgeo| = "
        f"{diff:.6e} m"
    )


    if diff > 1.0e-6:

        initial_ok = False


print()


# ------------------------------------------------------------
# TARGET-TIME SUMMARY
# ------------------------------------------------------------

summary = {}


print("=" * 82)
print("TIME-SLICE SUMMARY")
print("=" * 82)


for key in CASE_ORDER:

    res = results[key]

    summary[key] = {
        "time": [],
        "dzmax": [],
        "eta": [],
        "q": [],
        "max_raise": [],
        "min_bottom": [],
    }


    print()
    print(
        f"{key} | "
        f"k_d = {CASES[key]['kd']:.1e} m/(Pa s)"
    )

    print("-" * 82)


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


        max_raise = float(
            np.max(
                z - z0_geo
            )
        )


        min_bottom = float(
            np.min(z)
        )


        summary[key]["time"].append(
            actual_time
        )

        summary[key]["dzmax"].append(
            dzmax
        )

        summary[key]["eta"].append(
            eta
        )

        summary[key]["q"].append(
            q
        )

        summary[key]["max_raise"].append(
            max_raise
        )

        summary[key]["min_bottom"].append(
            min_bottom
        )


        print(
            f"t={actual_time:6.2f} s | "
            f"max lowering={dzmax:9.6f} m | "
            f"eta_res={eta:9.6f} m | "
            f"Qcrest={q:11.8f} m3/s"
        )


    print()

    print(
        "Maximum raising diagnostic = "
        f"{max(summary[key]['max_raise']):.6e} m"
    )

    print(
        "Minimum BOTTOM = "
        f"{min(summary[key]['min_bottom']):.6f} m"
    )


# ------------------------------------------------------------
# FULL TIME SERIES
# ------------------------------------------------------------

print()
print("=" * 82)
print("BUILDING FULL TIME SERIES")
print("=" * 82)


full_series = {}


for key in CASE_ORDER:

    res = results[key]

    times = np.asarray(
        res.times,
        dtype=float
    )


    dz_series = []
    eta_series = []
    q_series = []


    for rec in range(
        len(times)
    ):

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


    full_series[key] = {
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
        f"{key}: "
        f"{len(times)} records processed."
    )


# ------------------------------------------------------------
# PHYSICAL / NUMERICAL CONSISTENCY CHECK
# ------------------------------------------------------------

print()
print("=" * 82)
print("PHYSICAL / NUMERICAL CONSISTENCY CHECK")
print("=" * 82)


overall_pass = True


for key in CASE_ORDER:

    res = results[key]


    max_raise_all = -np.inf
    min_bottom_all = np.inf

    previous_z = None
    temporal_raise = 0.0


    for rec in range(
        len(res.times)
    ):

        z = get_bottom(
            res,
            rec
        )


        max_raise_all = max(
            max_raise_all,
            float(
                np.max(
                    z - z0_geo
                )
            )
        )


        min_bottom_all = min(
            min_bottom_all,
            float(
                np.min(z)
            )
        )


        if previous_z is not None:

            temporal_raise = max(
                temporal_raise,
                float(
                    np.max(
                        z - previous_z
                    )
                )
            )


        previous_z = z.copy()


    no_raising = (
        max_raise_all
        <= 1.0e-6
    )


    no_temporal_raising = (
        temporal_raise
        <= 1.0e-6
    )


    no_below_zero = (
        min_bottom_all
        >= -1.0e-6
    )


    case_pass = (
        no_raising
        and no_temporal_raising
        and no_below_zero
    )


    if not case_pass:

        overall_pass = False


    print(
        f"{key:6s} | "
        f"max raise={max_raise_all:.3e} m | "
        f"temporal raise={temporal_raise:.3e} m | "
        f"min Z={min_bottom_all:.6f} m | "
        f"{'PASS' if case_pass else 'FAIL'}"
    )


# ------------------------------------------------------------
# 30-S PARAMETER-RESPONSE TABLE
# ------------------------------------------------------------

print()
print("=" * 82)
print("30-S K_D RESPONSE")
print("=" * 82)

print(
    f"{'Case':<8}"
    f"{'k_d':>14}"
    f"{'Max erosion':>16}"
    f"{'Eta_res':>14}"
    f"{'Qcrest':>16}"
)

print("-" * 82)


for key in CASE_ORDER:

    kd = CASES[key]["kd"]

    dz30 = summary[key]["dzmax"][-1]
    eta30 = summary[key]["eta"][-1]
    q30 = summary[key]["q"][-1]


    print(
        f"{key:<8}"
        f"{kd:>14.3e}"
        f"{dz30:>16.6f}"
        f"{eta30:>14.6f}"
        f"{q30:>16.8f}"
    )


# ------------------------------------------------------------
# FIGURE 1
# MAXIMUM EROSION VS TIME
# ------------------------------------------------------------

plt.figure(
    figsize=(8.5, 5.8)
)


for key in CASE_ORDER:

    data = full_series[key]

    plt.plot(
        data["time"],
        data["dzmax"],
        linewidth=2.0,
        label=CASES[key]["label"]
    )


plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    "Maximum bed lowering (m)"
)

plt.title(
    "B-2 sensitivity of maximum erosion to erodibility coefficient"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()


FIG1 = os.path.join(
    FIG_DIR,
    "B2_kd_max_erosion_vs_time.png"
)


plt.savefig(
    FIG1,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FIGURE 2
# FINAL LONGITUDINAL BED PROFILE
# ------------------------------------------------------------

plt.figure(
    figsize=(9.5, 5.8)
)


ref_res = results["KD500"]

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


for key in CASE_ORDER:

    res = results[key]

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
        label=CASES[key]["label"]
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


FIG2 = os.path.join(
    FIG_DIR,
    "B2_kd_final_bed_profiles.png"
)


plt.savefig(
    FIG2,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FIGURE 3
# RESERVOIR WATER LEVEL
# ------------------------------------------------------------

plt.figure(
    figsize=(8.5, 5.8)
)


for key in CASE_ORDER:

    data = full_series[key]

    plt.plot(
        data["time"],
        data["eta"],
        linewidth=2.0,
        label=CASES[key]["label"]
    )


plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    "Reservoir water level (m)"
)

plt.title(
    "B-2 reservoir-level response to erodibility coefficient"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()


FIG3 = os.path.join(
    FIG_DIR,
    "B2_kd_reservoir_level.png"
)


plt.savefig(
    FIG3,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FIGURE 4
# CREST DISCHARGE
# ------------------------------------------------------------

plt.figure(
    figsize=(8.5, 5.8)
)


for key in CASE_ORDER:

    data = full_series[key]

    plt.plot(
        data["time"],
        data["q"],
        linewidth=2.0,
        label=CASES[key]["label"]
    )


plt.xlabel(
    "Time (s)"
)

plt.ylabel(
    r"Crest discharge $Q$ ($m^3/s$)"
)

plt.title(
    "B-2 crest-discharge response to erodibility coefficient"
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()


FIG4 = os.path.join(
    FIG_DIR,
    "B2_kd_crest_discharge.png"
)


plt.savefig(
    FIG4,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FIGURE 5
# PARAMETER RESPONSE AT 30 S
# ------------------------------------------------------------

kd_values = np.asarray(
    [
        CASES[key]["kd"]
        for key in CASE_ORDER
    ],
    dtype=float
)


dz30_values = np.asarray(
    [
        summary[key]["dzmax"][-1]
        for key in CASE_ORDER
    ],
    dtype=float
)


plt.figure(
    figsize=(7.5, 5.5)
)


plt.plot(
    kd_values,
    dz30_values,
    marker="o",
    linewidth=2.0
)


plt.xlabel(
    r"$k_d$ (m/(Pa s))"
)

plt.ylabel(
    "Maximum bed lowering at 30 s (m)"
)

plt.title(
    "B-2 erosion response to erodibility coefficient"
)

plt.grid(
    True,
    alpha=0.3
)

plt.tight_layout()


FIG5 = os.path.join(
    FIG_DIR,
    "B2_kd_parameter_response_30s.png"
)


plt.savefig(
    FIG5,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# FINAL REPORT
# ------------------------------------------------------------

print()
print("=" * 82)
print("FIGURES SAVED")
print("=" * 82)


for path in [
    FIG1,
    FIG2,
    FIG3,
    FIG4,
    FIG5,
]:

    print(
        os.path.basename(path)
    )


print()
print("=" * 82)


if (
    initial_ok
    and overall_pass
):

    print(
        "STATUS: B2 FIVE-CASE K_D ANALYSIS PASSED"
    )

else:

    print(
        "STATUS: B2 FIVE-CASE K_D ANALYSIS REQUIRES REVIEW"
    )


print("=" * 82)


# ------------------------------------------------------------
# CLOSE FILES
# ------------------------------------------------------------

for res in results.values():

    res.close()


geo.close()
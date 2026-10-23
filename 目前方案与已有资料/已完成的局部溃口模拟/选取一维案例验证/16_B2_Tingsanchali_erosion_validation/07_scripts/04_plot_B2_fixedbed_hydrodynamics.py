# -*- coding: utf-8 -*-

"""
B2 - Tingsanchali & Chinnarasri (2001), Run B-2

Generate fixed-bed hydrodynamic result figures directly from TELEMAC SLF.

Outputs:
1. B2_fixedbed_longitudinal_profile.png
2. B2_fixedbed_reservoir_level.png
3. B2_fixedbed_crest_discharge.png

Physical coordinate:
    x = 0.00 m   : upstream/front edge of dam crest
    x = 0.30 m   : downstream/back edge of dam crest
    x = -2.40 m  : upstream toe
    x = 4.30 m   : downstream toe

The longitudinal profile is extracted along the flume centreline:
    y = 0.50 m
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from data_manip.extraction.telemac_file import TelemacFile


# ============================================================
# 1. PATHS
# ============================================================

ROOT = (
    r"F:\Telemac_Projects\Tuanzho_New"
    r"\16_B2_Tingsanchali_erosion_validation"
)

RES_FILE = os.path.join(
    ROOT,
    "05_run",
    "B2_Tingsanchali_fixedbed_res.slf",
)

FIG_DIR = os.path.join(
    ROOT,
    "08_figures",
)

os.makedirs(
    FIG_DIR,
    exist_ok=True,
)


# ============================================================
# 2. OUTPUT FILES
# ============================================================

PROFILE_FIG = os.path.join(
    FIG_DIR,
    "B2_fixedbed_longitudinal_profile.png",
)

LEVEL_FIG = os.path.join(
    FIG_DIR,
    "B2_fixedbed_reservoir_level.png",
)

DISCHARGE_FIG = os.path.join(
    FIG_DIR,
    "B2_fixedbed_crest_discharge.png",
)


# ============================================================
# 3. PHYSICAL / EXTRACTION SETTINGS
# ============================================================

# Full flume
X_MIN = -17.50
X_MAX = 17.50

# Centreline
Y_CENTER = 0.50

# Dam geometry
X_UPSTREAM_TOE = -2.40
X_CREST_FRONT = 0.00
X_CREST_BACK = 0.30
X_DOWNSTREAM_TOE = 4.30

# Initial reservoir level
ETA_INITIAL = 0.85

# Representative upstream reservoir region
X_RES_MIN = -10.0
X_RES_MAX = -5.0

# Times shown in longitudinal profile
PROFILE_TIMES = [
    0.0,
    5.0,
    10.0,
    20.0,
    30.0,
]

# Longitudinal profile sampling
N_PROFILE_SEGMENTS = 700

# Crest cross-section sampling
N_CREST_SEGMENTS = 100

# Dry tolerance for plotting water surface
H_DRY = 1.0e-4


# ============================================================
# 4. OPEN TELEMAC RESULT
# ============================================================

if not os.path.isfile(RES_FILE):

    raise FileNotFoundError(
        f"Result file not found:\n{RES_FILE}"
    )


res = TelemacFile(RES_FILE)


try:

    times = np.asarray(
        res.times,
        dtype=float,
    )

    meshx = np.asarray(
        res.meshx,
        dtype=float,
    )

    meshy = np.asarray(
        res.meshy,
        dtype=float,
    )


    print("")
    print("=" * 78)
    print("B2 FIXED-BED HYDRODYNAMIC FIGURE GENERATOR")
    print("=" * 78)
    print(f"Result file       : {RES_FILE}")
    print(f"Records           : {len(times)}")
    print(f"Time range        : {times[0]:.2f} - {times[-1]:.2f} s")


    # ========================================================
    # 5. VARIABLE DETECTION
    # ========================================================

    variable_names = [
        str(v).strip()
        for v in res.varnames
    ]


    def find_variable(candidates):

        for candidate in candidates:

            for name in variable_names:

                if candidate.upper() in name.upper():

                    return name

        raise RuntimeError(
            "Could not find variable from candidates: "
            + ", ".join(candidates)
        )


    H_NAME = find_variable(
        [
            "WATER DEPTH",
            "DEPTH",
            "HAUTEUR",
        ]
    )

    U_NAME = find_variable(
        [
            "VELOCITY U",
            "U VELOCITY",
            "VITESSE U",
        ]
    )

    Z_NAME = find_variable(
        [
            "BOTTOM",
            "FOND",
        ]
    )


    print("")
    print("Detected variables:")
    print(f"  H                : {H_NAME}")
    print(f"  U                : {U_NAME}")
    print(f"  BOTTOM           : {Z_NAME}")


    # ========================================================
    # 6. HELPER FUNCTIONS
    # ========================================================

    def nearest_record(target_time):

        return int(
            np.argmin(
                np.abs(
                    times - target_time
                )
            )
        )


    def extract_polyline(variable_name, p1, p2, n_segments):

        points = [
            list(p1),
            list(p2),
        ]

        polyline_num = [
            n_segments,
        ]

        _, abs_curv, values = (
            res.get_timeseries_on_polyline(
                variable_name,
                points,
                polyline_num,
            )
        )

        return (
            np.asarray(abs_curv, dtype=float),
            np.asarray(values, dtype=float),
        )


    # ========================================================
    # 7. LONGITUDINAL CENTRELINE EXTRACTION
    # ========================================================

    print("")
    print("Extracting longitudinal centreline...")

    profile_p1 = [
        X_MIN,
        Y_CENTER,
    ]

    profile_p2 = [
        X_MAX,
        Y_CENTER,
    ]


    s_profile, H_profile = extract_polyline(
        H_NAME,
        profile_p1,
        profile_p2,
        N_PROFILE_SEGMENTS,
    )

    _, Z_profile = extract_polyline(
        Z_NAME,
        profile_p1,
        profile_p2,
        N_PROFILE_SEGMENTS,
    )


    # Since the line is horizontal and begins at X_MIN:
    x_profile = (
        X_MIN
        + s_profile
    )


    # ========================================================
    # 8. FIGURE 1:
    #    LONGITUDINAL WATER-SURFACE PROFILE
    # ========================================================

    print("Creating longitudinal profile figure...")


    fig, ax = plt.subplots(
        figsize=(11.0, 5.8)
    )


    # Bottom is fixed, so use first record
    bottom = Z_profile[:, 0]


    ax.plot(
        x_profile,
        bottom,
        linewidth=2.2,
        label="Bed elevation",
    )


    for target_time in PROFILE_TIMES:

        record = nearest_record(
            target_time
        )

        actual_time = times[record]

        h_now = H_profile[:, record]

        z_now = Z_profile[:, record]

        eta_now = (
            h_now
            + z_now
        )


        # Do not draw a false "water surface" over dry cells.
        eta_plot = eta_now.copy()

        eta_plot[
            h_now <= H_DRY
        ] = np.nan


        ax.plot(
            x_profile,
            eta_plot,
            linewidth=1.7,
            label=f"Water surface, t = {actual_time:.0f} s",
        )


    # Mark key dam positions
    ax.axvline(
        X_UPSTREAM_TOE,
        linewidth=0.8,
        linestyle="--",
    )

    ax.axvline(
        X_CREST_FRONT,
        linewidth=0.8,
        linestyle="--",
    )

    ax.axvline(
        X_CREST_BACK,
        linewidth=0.8,
        linestyle="--",
    )

    ax.axvline(
        X_DOWNSTREAM_TOE,
        linewidth=0.8,
        linestyle="--",
    )


    ax.set_xlabel(
        "Longitudinal coordinate x (m)"
    )

    ax.set_ylabel(
        "Elevation (m)"
    )

    ax.set_title(
        "B-2 fixed-bed longitudinal water-surface profiles"
    )

    ax.set_xlim(
        X_MIN,
        X_MAX,
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend(
        loc="best",
        fontsize=8.5,
    )

    fig.tight_layout()

    fig.savefig(
        PROFILE_FIG,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


    # ========================================================
    # 9. UPSTREAM RESERVOIR LEVEL TIME SERIES
    # ========================================================

    print("Calculating upstream reservoir level...")


    reservoir_mask = (
        (meshx >= X_RES_MIN - 1.0e-10)
        & (meshx <= X_RES_MAX + 1.0e-10)
    )


    if np.count_nonzero(reservoir_mask) == 0:

        raise RuntimeError(
            "No nodes found in reservoir sampling region."
        )


    eta_res_series = np.zeros(
        len(times),
        dtype=float,
    )


    for record in range(len(times)):

        H = np.asarray(
            res.get_data_value(
                H_NAME,
                record,
            ),
            dtype=float,
        )

        ZF = np.asarray(
            res.get_data_value(
                Z_NAME,
                record,
            ),
            dtype=float,
        )

        eta = H + ZF

        eta_res_series[record] = np.mean(
            eta[reservoir_mask]
        )


    # ========================================================
    # 10. FIGURE 2:
    #     RESERVOIR LEVEL
    # ========================================================

    print("Creating reservoir-level figure...")


    fig, ax = plt.subplots(
        figsize=(8.0, 5.0)
    )


    ax.plot(
        times,
        eta_res_series,
        linewidth=2.0,
        marker="o",
        markersize=3.5,
        label="Computed upstream water level",
    )


    ax.axhline(
        ETA_INITIAL,
        linewidth=1.0,
        linestyle="--",
        label="Initial water level = 0.85 m",
    )


    ax.set_xlabel(
        "Time (s)"
    )

    ax.set_ylabel(
        "Upstream water level (m)"
    )

    ax.set_title(
        "B-2 fixed-bed upstream reservoir level"
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend(
        loc="best"
    )

    fig.tight_layout()

    fig.savefig(
        LEVEL_FIG,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


    # ========================================================
    # 11. CREST DISCHARGE TIME SERIES
    # ========================================================

    print("Calculating crest discharge...")


    crest_p1 = [
        X_CREST_FRONT,
        0.0,
    ]

    crest_p2 = [
        X_CREST_FRONT,
        1.0,
    ]


    s_crest, H_crest = extract_polyline(
        H_NAME,
        crest_p1,
        crest_p2,
        N_CREST_SEGMENTS,
    )

    _, U_crest = extract_polyline(
        U_NAME,
        crest_p1,
        crest_p2,
        N_CREST_SEGMENTS,
    )


    q_crest_series = np.zeros(
        len(times),
        dtype=float,
    )


    for record in range(len(times)):

        h_now = H_crest[:, record]

        u_now = U_crest[:, record]

        unit_discharge = (
            h_now
            * u_now
        )

        q_crest_series[record] = np.trapz(
            unit_discharge,
            s_crest,
        )


    # ========================================================
    # 12. FIGURE 3:
    #     CREST DISCHARGE
    # ========================================================

    print("Creating crest-discharge figure...")


    fig, ax = plt.subplots(
        figsize=(8.0, 5.0)
    )


    ax.plot(
        times,
        q_crest_series,
        linewidth=2.0,
        marker="o",
        markersize=3.5,
        label="Computed discharge at x = 0 m",
    )


    ax.set_xlabel(
        "Time (s)"
    )

    ax.set_ylabel(
        r"Discharge $Q$ ($m^3/s$)"
    )

    ax.set_title(
        "B-2 fixed-bed overtopping discharge"
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend(
        loc="best"
    )

    fig.tight_layout()

    fig.savefig(
        DISCHARGE_FIG,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


    # ========================================================
    # 13. NUMERICAL SUMMARY
    # ========================================================

    q_peak_index = int(
        np.argmax(
            q_crest_series
        )
    )


    print("")
    print("=" * 78)
    print("HYDRODYNAMIC FIGURE SUMMARY")
    print("=" * 78)

    print(
        f"Initial reservoir level : "
        f"{eta_res_series[0]:.6f} m"
    )

    print(
        f"Final reservoir level   : "
        f"{eta_res_series[-1]:.6f} m"
    )

    print(
        f"Peak crest discharge    : "
        f"{q_crest_series[q_peak_index]:.8f} m3/s"
    )

    print(
        f"Peak-Q time             : "
        f"{times[q_peak_index]:.2f} s"
    )

    print(
        f"Final crest discharge   : "
        f"{q_crest_series[-1]:.8f} m3/s"
    )


    print("")
    print("Figures written:")

    print(
        f"  {PROFILE_FIG}"
    )

    print(
        f"  {LEVEL_FIG}"
    )

    print(
        f"  {DISCHARGE_FIG}"
    )

    print("")
    print(
        "STATUS: B2 FIXED-BED HYDRODYNAMIC "
        "FIGURES GENERATED"
    )

    print("=" * 78)


finally:

    res.close()
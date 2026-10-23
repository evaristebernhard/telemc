# -*- coding: utf-8 -*-

"""
B2 - Tingsanchali & Chinnarasri (2001), Run B-2
Check corrected full-flume fixed-bed hydrodynamics.

Checks:
1. SLF time records
2. NaN / Inf
3. Negative water depth
4. Maximum water depth and velocity
5. Upstream reservoir water level
6. Discharge across the dam crest section x = 0
7. Downstream wet-front position

This is a FIXED-BED hydrodynamic check.
No erosion validation is performed here.
"""

import os
import sys
import numpy as np

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


# ============================================================
# 2. PHYSICAL PARAMETERS
# ============================================================

ETA0 = 0.85

X_CREST_FRONT = 0.0

# Reservoir-level sampling region.
# Stay sufficiently upstream of the dam to avoid the strong local
# drawdown immediately around the crest.
X_RES_MIN = -10.0
X_RES_MAX = -5.0

# Numerical wet/dry threshold used only for diagnostics.
H_WET = 1.0e-4

# Cross-section for overtopping discharge.
# The mesh spacing is 0.05 m, so x = 0 is represented directly.
SECTION_P1 = [0.0, 0.0]
SECTION_P2 = [0.0, 1.0]

CHECK_TIMES = [
    0.0,
    5.0,
    10.0,
    20.0,
    30.0,
]


# ============================================================
# 3. OPEN RESULT FILE
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


    # ========================================================
    # 4. BASIC FILE INFORMATION
    # ========================================================

    print("")
    print("=" * 78)
    print(
        "B2 TINGSANCHALI FULL-FLUME FIXED-BED "
        "HYDRODYNAMIC CHECK"
    )
    print("=" * 78)

    print(
        f"Result file            : {RES_FILE}"
    )

    print(
        f"Mesh nodes             : {meshx.size}"
    )

    print(
        f"Number of records      : {len(times)}"
    )

    print(
        f"First time             : {times[0]:.6f} s"
    )

    print(
        f"Last time              : {times[-1]:.6f} s"
    )


    if meshx.size != 14721:

        raise RuntimeError(
            f"Unexpected node count: {meshx.size}; "
            "expected 14721."
        )


    if abs(times[0]) > 1.0e-10:

        raise RuntimeError(
            "First result time is not 0 s."
        )


    if abs(times[-1] - 30.0) > 1.0e-6:

        raise RuntimeError(
            "Last result time is not 30 s."
        )


    # ========================================================
    # 5. VARIABLE NAMES
    # ========================================================

    variable_names = [
        str(v).strip()
        for v in res.varnames
    ]

    print("")
    print("Variables:")
    for name in variable_names:
        print(f"  {name}")


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

    V_NAME = find_variable(
        [
            "VELOCITY V",
            "V VELOCITY",
            "VITESSE V",
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
    print(f"  H                     : {H_NAME}")
    print(f"  U                     : {U_NAME}")
    print(f"  V                     : {V_NAME}")
    print(f"  BOTTOM                : {Z_NAME}")


    # ========================================================
    # 6. RESERVOIR SAMPLING MASK
    # ========================================================

    reservoir_mask = (
        (meshx >= X_RES_MIN - 1.0e-10)
        & (meshx <= X_RES_MAX + 1.0e-10)
    )


    if np.count_nonzero(reservoir_mask) == 0:

        raise RuntimeError(
            "No nodes found in reservoir sampling region."
        )


    print("")
    print(
        "Reservoir sampling region:"
    )

    print(
        f"  x = {X_RES_MIN:.2f} to "
        f"{X_RES_MAX:.2f} m"
    )

    print(
        f"  nodes                 : "
        f"{np.count_nonzero(reservoir_mask)}"
    )


    # ========================================================
    # 7. HELPER: NEAREST RECORD
    # ========================================================

    def nearest_record(target_time):

        idx = int(
            np.argmin(
                np.abs(times - target_time)
            )
        )

        return idx


    # ========================================================
    # 8. HELPER: DISCHARGE THROUGH x = 0
    # ========================================================

    def get_crest_discharge(record):

        """
        Use TELEMAC extraction along the vertical section x = 0.

        The returned polyline is sampled densely enough for
        numerical integration of q = H * U across the 1 m width.
        """

        npoints = 101

        polyline_points = [
            SECTION_P1,
            SECTION_P2,
        ]

        polyline_num = [
            npoints - 1,
        ]

        _, abs_curv, h_section = (
            res.get_timeseries_on_polyline(
                H_NAME,
                polyline_points,
                polyline_num,
            )
        )

        _, _, u_section = (
            res.get_timeseries_on_polyline(
                U_NAME,
                polyline_points,
                polyline_num,
            )
        )

        h_now = np.asarray(
            h_section[:, record],
            dtype=float,
        )

        u_now = np.asarray(
            u_section[:, record],
            dtype=float,
        )

        s = np.asarray(
            abs_curv,
            dtype=float,
        )

        q_unit = h_now * u_now

        discharge = np.trapz(
            q_unit,
            s,
        )

        return float(discharge)


    # ========================================================
    # 9. CHECK SELECTED TIMES
    # ========================================================

    print("")
    print("-" * 78)

    print(
        "Selected-time hydrodynamic diagnostics"
    )

    print("-" * 78)

    print(
        f"{'t(s)':>7s} "
        f"{'eta_res(m)':>12s} "
        f"{'Qcrest(m3/s)':>14s} "
        f"{'maxH(m)':>11s} "
        f"{'maxV(m/s)':>12s} "
        f"{'wet_front_x(m)':>15s}"
    )


    global_min_h = np.inf
    global_max_h = -np.inf
    global_max_v = -np.inf

    found_nan = False
    found_inf = False

    selected_results = []


    for target_time in CHECK_TIMES:

        record = nearest_record(
            target_time
        )

        actual_time = times[record]


        H = np.asarray(
            res.get_data_value(
                H_NAME,
                record,
            ),
            dtype=float,
        )

        U = np.asarray(
            res.get_data_value(
                U_NAME,
                record,
            ),
            dtype=float,
        )

        V = np.asarray(
            res.get_data_value(
                V_NAME,
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


        if (
            np.isnan(H).any()
            or np.isnan(U).any()
            or np.isnan(V).any()
            or np.isnan(ZF).any()
        ):

            found_nan = True


        if (
            np.isinf(H).any()
            or np.isinf(U).any()
            or np.isinf(V).any()
            or np.isinf(ZF).any()
        ):

            found_inf = True


        velocity = np.sqrt(
            U ** 2
            + V ** 2
        )


        eta = H + ZF


        eta_res = float(
            np.mean(
                eta[reservoir_mask]
            )
        )


        max_h = float(
            np.nanmax(H)
        )

        min_h = float(
            np.nanmin(H)
        )

        max_v = float(
            np.nanmax(velocity)
        )


        global_min_h = min(
            global_min_h,
            min_h,
        )

        global_max_h = max(
            global_max_h,
            max_h,
        )

        global_max_v = max(
            global_max_v,
            max_v,
        )


        wet_mask = (
            (H > H_WET)
            & (meshx > X_CREST_FRONT)
        )


        if np.any(wet_mask):

            wet_front_x = float(
                np.max(
                    meshx[wet_mask]
                )
            )

        else:

            wet_front_x = np.nan


        q_crest = get_crest_discharge(
            record
        )


        selected_results.append(
            (
                actual_time,
                eta_res,
                q_crest,
                max_h,
                max_v,
                wet_front_x,
            )
        )


        wet_front_text = (
            f"{wet_front_x:15.6f}"
            if np.isfinite(wet_front_x)
            else f"{'DRY':>15s}"
        )


        print(
            f"{actual_time:7.2f} "
            f"{eta_res:12.6f} "
            f"{q_crest:14.8f} "
            f"{max_h:11.6f} "
            f"{max_v:12.6f} "
            f"{wet_front_text}"
        )


    # ========================================================
    # 10. SCAN ALL RECORDS FOR NUMERICAL PROBLEMS
    # ========================================================

    print("")
    print("-" * 78)
    print(
        "Scanning all records for numerical problems..."
    )
    print("-" * 78)


    all_min_h = np.inf
    all_max_h = -np.inf
    all_max_v = -np.inf

    all_nan = False
    all_inf = False


    for record in range(len(times)):

        H = np.asarray(
            res.get_data_value(
                H_NAME,
                record,
            ),
            dtype=float,
        )

        U = np.asarray(
            res.get_data_value(
                U_NAME,
                record,
            ),
            dtype=float,
        )

        V = np.asarray(
            res.get_data_value(
                V_NAME,
                record,
            ),
            dtype=float,
        )


        if (
            np.isnan(H).any()
            or np.isnan(U).any()
            or np.isnan(V).any()
        ):

            all_nan = True


        if (
            np.isinf(H).any()
            or np.isinf(U).any()
            or np.isinf(V).any()
        ):

            all_inf = True


        velocity = np.sqrt(
            U ** 2
            + V ** 2
        )


        all_min_h = min(
            all_min_h,
            float(np.nanmin(H)),
        )

        all_max_h = max(
            all_max_h,
            float(np.nanmax(H)),
        )

        all_max_v = max(
            all_max_v,
            float(np.nanmax(velocity)),
        )


    # ========================================================
    # 11. INITIAL CONDITION CHECK
    # ========================================================

    record0 = nearest_record(0.0)

    H0 = np.asarray(
        res.get_data_value(
            H_NAME,
            record0,
        ),
        dtype=float,
    )

    ZF0 = np.asarray(
        res.get_data_value(
            Z_NAME,
            record0,
        ),
        dtype=float,
    )

    ETA0_FIELD = H0 + ZF0


    upstream_initial_mask = (
        meshx <= 0.0 + 1.0e-10
    )

    upstream_wet_mask = (
        upstream_initial_mask
        & (H0 > H_WET)
    )


    if np.any(upstream_wet_mask):

        initial_eta_mean = float(
            np.mean(
                ETA0_FIELD[
                    upstream_wet_mask
                ]
            )
        )

        initial_eta_min = float(
            np.min(
                ETA0_FIELD[
                    upstream_wet_mask
                ]
            )
        )

        initial_eta_max = float(
            np.max(
                ETA0_FIELD[
                    upstream_wet_mask
                ]
            )
        )

    else:

        raise RuntimeError(
            "No initially wet upstream nodes found."
        )


    downstream_initial_wet = np.count_nonzero(
        (meshx > 0.0 + 1.0e-10)
        & (H0 > H_WET)
    )


    # ========================================================
    # 12. FINAL REPORT
    # ========================================================

    print("")
    print("=" * 78)

    print(
        "FIXED-BED HYDRODYNAMIC SUMMARY"
    )

    print("=" * 78)


    print("")
    print("Initial condition:")

    print(
        f"  Mean upstream eta     : "
        f"{initial_eta_mean:.8f} m"
    )

    print(
        f"  Min upstream eta      : "
        f"{initial_eta_min:.8f} m"
    )

    print(
        f"  Max upstream eta      : "
        f"{initial_eta_max:.8f} m"
    )

    print(
        f"  Initially wet nodes "
        f"downstream of x=0 : "
        f"{downstream_initial_wet}"
    )


    print("")
    print("Whole-run numerical range:")

    print(
        f"  Minimum H             : "
        f"{all_min_h:.12e} m"
    )

    print(
        f"  Maximum H             : "
        f"{all_max_h:.8f} m"
    )

    print(
        f"  Maximum velocity      : "
        f"{all_max_v:.8f} m/s"
    )

    print(
        f"  NaN detected          : "
        f"{all_nan}"
    )

    print(
        f"  Inf detected          : "
        f"{all_inf}"
    )


    # ========================================================
    # 13. PASS / WARNING LOGIC
    # ========================================================

    passed = True


    if all_nan:

        print("")
        print(
            "FAIL: NaN values detected."
        )

        passed = False


    if all_inf:

        print("")
        print(
            "FAIL: Inf values detected."
        )

        passed = False


    if all_min_h < -1.0e-6:

        print("")
        print(
            "FAIL: Significant negative water "
            "depth detected."
        )

        passed = False


    if abs(
        initial_eta_mean
        - ETA0
    ) > 1.0e-6:

        print("")
        print(
            "FAIL: Initial upstream water level "
            "does not match 0.85 m."
        )

        passed = False


    if downstream_initial_wet != 0:

        print("")
        print(
            "WARNING: Some nodes downstream of "
            "x=0 are initially wet."
        )


    print("")
    print("=" * 78)


    if passed:

        print(
            "STATUS: FIXED-BED HYDRODYNAMIC "
            "CHECK PASSED"
        )

    else:

        print(
            "STATUS: FIXED-BED HYDRODYNAMIC "
            "CHECK FAILED"
        )


    print("=" * 78)


finally:

    res.close()
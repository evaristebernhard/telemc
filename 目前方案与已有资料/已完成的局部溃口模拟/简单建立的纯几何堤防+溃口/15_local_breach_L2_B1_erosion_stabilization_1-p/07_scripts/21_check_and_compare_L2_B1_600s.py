# -*- coding: utf-8 -*-

"""
TUANZHOU L2-B1 600 s technical check and B0/B1 comparison

Checks:
1. B1 result reaches 600 s.
2. No bed raising.
3. No erosion outside the prescribed activation zone.
4. No BOTTOM below ZMIN.
5. Cumulative maximum erosion is non-decreasing.
6. Number of cumulatively eroded nodes is non-decreasing.
7. Quantitative B0/B1 comparison at:
   60, 120, 180, 300, 420, 600 s.

Cumulative erosion:
    erosion = BOTTOM(t=0) - BOTTOM(t)

Positive erosion means bed lowering.
"""

from pathlib import Path
import sys
import numpy as np


# ======================================================================
# PATHS
# ======================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

B1_FILE = (
    PROJECT_DIR
    / "05_run"
    / "local_breach_L2_B1_600s_res.slf"
)

B0_FILE = Path(
    r"F:\Telemac_Projects\Tuanzho_New"
    r"\14_local_breach_L2_B0_vertical_erosion"
    r"\09_B0_600s_test"
    r"\05_run"
    r"\local_breach_L2_B0_600s_res.slf"
)


# ======================================================================
# PARAMETERS
# ======================================================================

CHECK_TIMES = [
    60.0,
    120.0,
    180.0,
    300.0,
    420.0,
    600.0,
]

EROSION_TOL = 1.0e-6

ZMIN = 27.0

X_EROSION_MIN = -113.0
X_EROSION_MAX = 113.0

Y_EROSION_MIN = -25.0
Y_EROSION_MAX = 35.0

INITIAL_BREACH_LEFT = -10.0
INITIAL_BREACH_RIGHT = 10.0

MESH_TOL = 1.0e-6


# ======================================================================
# TELEMAC IMPORT
# ======================================================================

try:
    from data_manip.extraction.telemac_file import TelemacFile
except ImportError:
    print()
    print("ERROR:")
    print("Could not import TELEMAC Python module:")
    print("data_manip.extraction.telemac_file")
    print()
    print("Run this script inside the configured openTELEMAC environment.")
    print()
    sys.exit(1)


# ======================================================================
# FUNCTIONS
# ======================================================================

def nearest_record(times, target_time):
    return int(
        np.argmin(
            np.abs(times - target_time)
        )
    )


def get_bottom(res, record):

    candidates = [
        "BOTTOM",
        "BOTTOM          ",
        "FOND",
    ]

    available = [
        str(v).strip()
        for v in res.varnames
    ]

    for name in candidates:

        if name.strip() in available:

            return np.asarray(
                res.get_data_value(
                    name.strip(),
                    record,
                ),
                dtype=float,
            )

    print()
    print("ERROR: BOTTOM variable was not found.")
    print("Available variables:")

    for name in available:
        print("   ", repr(name))

    print()
    sys.exit(1)


def cumulative_erosion(bottom0, bottom):

    erosion = bottom0 - bottom

    erosion[
        np.abs(erosion) < EROSION_TOL
    ] = 0.0

    return erosion


def analyse_record(
    x,
    y,
    bottom0,
    bottom,
):

    erosion = cumulative_erosion(
        bottom0,
        bottom,
    )

    eroded_mask = erosion > EROSION_TOL

    raised_mask = erosion < -EROSION_TOL

    activation_mask = (
        (x >= X_EROSION_MIN)
        & (x <= X_EROSION_MAX)
        & (y >= Y_EROSION_MIN)
        & (y <= Y_EROSION_MAX)
    )

    outside_activation = (
        eroded_mask
        & (~activation_mask)
    )

    inside_initial = (
        eroded_mask
        & (x >= INITIAL_BREACH_LEFT)
        & (x <= INITIAL_BREACH_RIGHT)
    )

    outside_initial = (
        eroded_mask
        & (
            (x < INITIAL_BREACH_LEFT)
            | (x > INITIAL_BREACH_RIGHT)
        )
    )

    if np.any(eroded_mask):

        strongest = int(
            np.argmax(erosion)
        )

        max_erosion = float(
            erosion[strongest]
        )

        strongest_x = float(
            x[strongest]
        )

        strongest_y = float(
            y[strongest]
        )

    else:

        strongest = -1
        max_erosion = 0.0
        strongest_x = np.nan
        strongest_y = np.nan

    return {
        "erosion": erosion,
        "max_erosion": max_erosion,
        "eroded_nodes": int(
            np.count_nonzero(eroded_mask)
        ),
        "raised_nodes": int(
            np.count_nonzero(raised_mask)
        ),
        "outside_activation": int(
            np.count_nonzero(outside_activation)
        ),
        "inside_initial": int(
            np.count_nonzero(inside_initial)
        ),
        "outside_initial": int(
            np.count_nonzero(outside_initial)
        ),
        "below_zmin": int(
            np.count_nonzero(
                bottom < (ZMIN - EROSION_TOL)
            )
        ),
        "strongest_index": strongest,
        "strongest_x": strongest_x,
        "strongest_y": strongest_y,
        "sum_erosion": float(
            np.sum(
                np.maximum(
                    erosion,
                    0.0,
                )
            )
        ),
    }


def print_separator():
    print("-" * 110)


# ======================================================================
# MAIN
# ======================================================================

def main():

    print()
    print("=" * 110)
    print("TUANZHOU L2-B1 600 s TECHNICAL CHECK + B0/B1 COMPARISON")
    print("=" * 110)
    print()

    print("B0 file:")
    print(B0_FILE)
    print()

    print("B1 file:")
    print(B1_FILE)
    print()

    if not B0_FILE.exists():

        print("ERROR: B0 result file does not exist.")
        sys.exit(1)

    if not B1_FILE.exists():

        print("ERROR: B1 result file does not exist.")
        sys.exit(1)

    b0 = TelemacFile(
        str(B0_FILE)
    )

    b1 = TelemacFile(
        str(B1_FILE)
    )

    try:

        # --------------------------------------------------------------
        # Basic data
        # --------------------------------------------------------------

        x0 = np.asarray(
            b0.meshx,
            dtype=float,
        )

        y0 = np.asarray(
            b0.meshy,
            dtype=float,
        )

        x1 = np.asarray(
            b1.meshx,
            dtype=float,
        )

        y1 = np.asarray(
            b1.meshy,
            dtype=float,
        )

        t0 = np.asarray(
            b0.times,
            dtype=float,
        )

        t1 = np.asarray(
            b1.times,
            dtype=float,
        )

        print(f"B0 nodes   : {len(x0)}")
        print(f"B1 nodes   : {len(x1)}")

        print(
            f"B0 records : {len(t0)}"
            f" | {t0[0]:.3f} -> {t0[-1]:.3f} s"
        )

        print(
            f"B1 records : {len(t1)}"
            f" | {t1[0]:.3f} -> {t1[-1]:.3f} s"
        )

        print()

        # --------------------------------------------------------------
        # Mesh consistency
        # --------------------------------------------------------------

        same_node_count = (
            len(x0) == len(x1)
        )

        same_mesh = False

        if same_node_count:

            max_dx = float(
                np.max(
                    np.abs(x0 - x1)
                )
            )

            max_dy = float(
                np.max(
                    np.abs(y0 - y1)
                )
            )

            same_mesh = (
                max_dx <= MESH_TOL
                and max_dy <= MESH_TOL
            )

        else:

            max_dx = np.nan
            max_dy = np.nan

        print("MESH CONSISTENCY")
        print_separator()

        print(
            f"Same node count : {same_node_count}"
        )

        print(
            f"Maximum |dx|    : {max_dx:.9e} m"
        )

        print(
            f"Maximum |dy|    : {max_dy:.9e} m"
        )

        print(
            f"Same mesh       : {same_mesh}"
        )

        print()

        if not same_mesh:

            print("ERROR:")
            print(
                "B0 and B1 meshes are not identical. "
                "Direct node-by-node comparison is invalid."
            )
            sys.exit(1)

        # --------------------------------------------------------------
        # Initial bottom consistency
        # --------------------------------------------------------------

        bottom0_b0 = get_bottom(
            b0,
            0,
        )

        bottom0_b1 = get_bottom(
            b1,
            0,
        )

        initial_bottom_difference = float(
            np.max(
                np.abs(
                    bottom0_b0
                    - bottom0_b1
                )
            )
        )

        print("INITIAL BOTTOM CONSISTENCY")
        print_separator()

        print(
            "Maximum |BOTTOM_B1(0)-BOTTOM_B0(0)| = "
            f"{initial_bottom_difference:.9e} m"
        )

        print()

        # --------------------------------------------------------------
        # B1 full-record technical checks
        # --------------------------------------------------------------

        max_series = []
        node_series = []

        total_raised_nodes = 0
        total_outside_activation = 0
        total_below_zmin = 0

        for rec in range(len(t1)):

            bottom = get_bottom(
                b1,
                rec,
            )

            stats = analyse_record(
                x1,
                y1,
                bottom0_b1,
                bottom,
            )

            max_series.append(
                stats["max_erosion"]
            )

            node_series.append(
                stats["eroded_nodes"]
            )

            total_raised_nodes = max(
                total_raised_nodes,
                stats["raised_nodes"],
            )

            total_outside_activation = max(
                total_outside_activation,
                stats["outside_activation"],
            )

            total_below_zmin = max(
                total_below_zmin,
                stats["below_zmin"],
            )

        max_series = np.asarray(
            max_series,
            dtype=float,
        )

        node_series = np.asarray(
            node_series,
            dtype=int,
        )

        max_monotonic = bool(
            np.all(
                np.diff(max_series)
                >= -EROSION_TOL
            )
        )

        node_monotonic = bool(
            np.all(
                np.diff(node_series)
                >= 0
            )
        )

        reaches_600 = bool(
            t1[-1] >= 599.999
        )

        no_raising = (
            total_raised_nodes == 0
        )

        no_outside_activation = (
            total_outside_activation == 0
        )

        no_below_zmin = (
            total_below_zmin == 0
        )

        print("B1 FULL-RECORD TECHNICAL CHECK")
        print_separator()

        print(
            f"Reaches 600 s                         : {reaches_600}"
        )

        print(
            f"No bed raising                        : {no_raising}"
        )

        print(
            f"No erosion outside activation zone    : "
            f"{no_outside_activation}"
        )

        print(
            f"No BOTTOM below ZMIN                  : {no_below_zmin}"
        )

        print(
            f"Maximum erosion monotonic             : {max_monotonic}"
        )

        print(
            f"Cumulative eroded-node count monotonic: {node_monotonic}"
        )

        print()

        # --------------------------------------------------------------
        # Selected-time comparison
        # --------------------------------------------------------------

        print("=" * 110)
        print("B0 / B1 SELECTED-TIME COMPARISON")
        print("=" * 110)

        header = (
            f"{'time':>7}"
            f"{'B0 max D':>14}"
            f"{'B1 max D':>14}"
            f"{'diff':>14}"
            f"{'rel.diff':>13}"
            f"{'B0 nodes':>11}"
            f"{'B1 nodes':>11}"
            f"{'B0 out20':>11}"
            f"{'B1 out20':>11}"
        )

        print(header)
        print_separator()

        comparison_rows = []

        for target_time in CHECK_TIMES:

            rec0 = nearest_record(
                t0,
                target_time,
            )

            rec1 = nearest_record(
                t1,
                target_time,
            )

            actual_t0 = float(
                t0[rec0]
            )

            actual_t1 = float(
                t1[rec1]
            )

            bottom_b0 = get_bottom(
                b0,
                rec0,
            )

            bottom_b1 = get_bottom(
                b1,
                rec1,
            )

            s0 = analyse_record(
                x0,
                y0,
                bottom0_b0,
                bottom_b0,
            )

            s1 = analyse_record(
                x1,
                y1,
                bottom0_b1,
                bottom_b1,
            )

            difference = (
                s1["max_erosion"]
                - s0["max_erosion"]
            )

            if s0["max_erosion"] > 0.0:

                relative = (
                    difference
                    / s0["max_erosion"]
                    * 100.0
                )

            else:

                relative = np.nan

            comparison_rows.append(
                {
                    "time": target_time,
                    "actual_t0": actual_t0,
                    "actual_t1": actual_t1,
                    "b0": s0,
                    "b1": s1,
                    "difference": difference,
                    "relative": relative,
                }
            )

            print(
                f"{target_time:7.0f}"
                f"{s0['max_erosion']:14.6f}"
                f"{s1['max_erosion']:14.6f}"
                f"{difference:14.6e}"
                f"{relative:12.3f}%"
                f"{s0['eroded_nodes']:11d}"
                f"{s1['eroded_nodes']:11d}"
                f"{s0['outside_initial']:11d}"
                f"{s1['outside_initial']:11d}"
            )

        print()

        # --------------------------------------------------------------
        # Strongest-point comparison
        # --------------------------------------------------------------

        print("=" * 110)
        print("STRONGEST EROSION POINT COMPARISON")
        print("=" * 110)

        for row in comparison_rows:

            target_time = row["time"]
            s0 = row["b0"]
            s1 = row["b1"]

            print()
            print(f"t = {target_time:.0f} s")

            print(
                "  B0: "
                f"Dmax={s0['max_erosion']:.9e} m"
                f" | x={s0['strongest_x']:.6f}"
                f" | y={s0['strongest_y']:.6f}"
            )

            print(
                "  B1: "
                f"Dmax={s1['max_erosion']:.9e} m"
                f" | x={s1['strongest_x']:.6f}"
                f" | y={s1['strongest_y']:.6f}"
            )

        print()

        # --------------------------------------------------------------
        # Final 600 s detailed check
        # --------------------------------------------------------------

        final_rec = nearest_record(
            t1,
            600.0,
        )

        final_bottom = get_bottom(
            b1,
            final_rec,
        )

        final_stats = analyse_record(
            x1,
            y1,
            bottom0_b1,
            final_bottom,
        )

        print("=" * 110)
        print("B1 FINAL 600 s SUMMARY")
        print("=" * 110)

        print(
            f"Actual time                  : "
            f"{t1[final_rec]:.3f} s"
        )

        print(
            f"Maximum cumulative erosion   : "
            f"{final_stats['max_erosion']:.9e} m"
        )

        print(
            f"Cumulatively eroded nodes    : "
            f"{final_stats['eroded_nodes']}"
        )

        print(
            f"Inside initial 20 m breach   : "
            f"{final_stats['inside_initial']}"
        )

        print(
            f"Outside initial 20 m breach  : "
            f"{final_stats['outside_initial']}"
        )

        print(
            f"Outside activation zone      : "
            f"{final_stats['outside_activation']}"
        )

        print(
            f"Raised nodes                 : "
            f"{final_stats['raised_nodes']}"
        )

        print(
            f"Nodes below ZMIN             : "
            f"{final_stats['below_zmin']}"
        )

        print(
            f"Strongest point x            : "
            f"{final_stats['strongest_x']:.6f} m"
        )

        print(
            f"Strongest point y            : "
            f"{final_stats['strongest_y']:.6f} m"
        )

        print(
            f"Sum of positive erosion      : "
            f"{final_stats['sum_erosion']:.9e} m"
        )

        print()

        # --------------------------------------------------------------
        # Final pass/fail
        # --------------------------------------------------------------

        technical_pass = all(
            [
                reaches_600,
                no_raising,
                no_outside_activation,
                no_below_zmin,
                max_monotonic,
                node_monotonic,
                same_mesh,
            ]
        )

        print("=" * 110)

        if technical_pass:

            print(
                "L2-B1 600 s LONG-DURATION TECHNICAL CHECK: PASS"
            )

        else:

            print(
                "L2-B1 600 s LONG-DURATION TECHNICAL CHECK: FAIL"
            )

        print("=" * 110)
        print()

        if not technical_pass:
            sys.exit(2)

    finally:

        b0.close()
        b1.close()


if __name__ == "__main__":
    main()
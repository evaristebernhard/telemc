import os
import sys
import numpy as np

# ============================================================
# B2 Tingsanchali erosion reference-run checker
#
# Purpose:
#   Validate the first porosity-inclusive excess-shear
#   erosion reference run before parameter calibration.
#
# Checks:
#   1. Result file and variables
#   2. NaN / Inf
#   3. Initial bottom consistency
#   4. Bed lowering with time
#   5. No bed raising
#   6. No erosion below Z = 0
#   7. No erosion outside initial dam footprint
#   8. Temporal accumulation of erosion
#
# Benchmark:
#   Tingsanchali & Chinnarasri (2001), Run B-2
#
# Reference erosion parameters:
#   p     = 0.40
#   tau_c = 0.01 Pa
#   k_d   = 5.0e-4 m/(Pa s)
#
# NOTE:
#   This script checks numerical / implementation behaviour.
#   It does NOT declare the erosion parameters calibrated.
# ============================================================


# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

ROOT = r"F:\Telemac_Projects\Tuanzho_New\16_B2_Tingsanchali_erosion_validation"

GEO_FILE = os.path.join(
    ROOT,
    "05_run",
    "B2_Tingsanchali_geo.slf"
)

RES_FILE = os.path.join(
    ROOT,
    "05_run",
    "B2_Tingsanchali_erosion_ref_res.slf"
)


# ------------------------------------------------------------
# IMPORT TELEMAC
# ------------------------------------------------------------

try:
    from data_manip.extraction.telemac_file import TelemacFile
except Exception as exc:
    print("ERROR: Could not import TelemacFile.")
    print(exc)
    sys.exit(1)


# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

CHECK_TIMES = [0.0, 5.0, 10.0, 20.0, 30.0]

TIME_TOL = 0.51

Z_TOL = 1.0e-8

CHANGE_TOL = 1.0e-7

X_DAM_MIN = -2.40
X_DAM_MAX = 4.30

Y_DAM_MIN = 0.00
Y_DAM_MAX = 1.00

Z_MIN = 0.0


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

def find_variable_name(varnames, candidates):

    stripped = [v.strip() for v in varnames]

    for candidate in candidates:
        for original, clean in zip(varnames, stripped):
            if clean.upper() == candidate.upper():
                return original

    return None


def nearest_record(times, target):

    idx = int(np.argmin(np.abs(times - target)))

    if abs(times[idx] - target) > TIME_TOL:
        raise RuntimeError(
            f"No output record close enough to t={target:.2f} s. "
            f"Nearest time is {times[idx]:.6f} s."
        )

    return idx


def print_separator():

    print("-" * 78)


# ------------------------------------------------------------
# FILE CHECK
# ------------------------------------------------------------

print()
print("=" * 78)
print("B2 TINGSANCHALI EROSION REFERENCE CHECK")
print("=" * 78)

if not os.path.isfile(GEO_FILE):
    print("ERROR: Geometry file not found:")
    print(GEO_FILE)
    sys.exit(1)

if not os.path.isfile(RES_FILE):
    print("ERROR: Result file not found:")
    print(RES_FILE)
    sys.exit(1)

print("Geometry :", GEO_FILE)
print("Result   :", RES_FILE)
print()


# ------------------------------------------------------------
# OPEN FILES
# ------------------------------------------------------------

geo = None
res = None

try:

    geo = TelemacFile(GEO_FILE)
    res = TelemacFile(RES_FILE)

    x = np.asarray(res.meshx, dtype=float)
    y = np.asarray(res.meshy, dtype=float)

    times = np.asarray(res.times, dtype=float)

    print(f"Mesh nodes       : {len(x)}")
    print(f"Result records   : {len(times)}")
    print(
        f"Result time span: "
        f"{times[0]:.6f} -> {times[-1]:.6f} s"
    )

    print()
    print("Result variables:")

    for name in res.varnames:
        print("  ", name.strip())

    print()


    # --------------------------------------------------------
    # FIND BOTTOM VARIABLE
    # --------------------------------------------------------

    bottom_geo_name = find_variable_name(
        geo.varnames,
        ["BOTTOM", "B"]
    )

    bottom_res_name = find_variable_name(
        res.varnames,
        ["BOTTOM", "B"]
    )

    if bottom_geo_name is None:
        raise RuntimeError(
            "BOTTOM variable not found in geometry file."
        )

    if bottom_res_name is None:
        raise RuntimeError(
            "BOTTOM variable not found in result file."
        )

    print(
        "Geometry bottom variable :",
        bottom_geo_name.strip()
    )

    print(
        "Result bottom variable   :",
        bottom_res_name.strip()
    )

    print()


    # --------------------------------------------------------
    # READ INITIAL GEOMETRY BOTTOM
    # --------------------------------------------------------

    z_initial_geometry = np.asarray(
        geo.get_data_value(bottom_geo_name, 0),
        dtype=float
    )


    # --------------------------------------------------------
    # DEFINE DAM / OUTSIDE-DAM MASKS
    #
    # IMPORTANT:
    #   The erosion Fortran uses the rectangular activation
    #   range -2.40 <= x <= 4.30 and 0 <= y <= 1.
    #
    #   However, the original flume floor inside this x-range
    #   already has Z = 0 at the toe endpoints and cannot
    #   erode because the Fortran additionally requires ZF>0.
    # --------------------------------------------------------

    dam_box_mask = (
        (x >= X_DAM_MIN - Z_TOL)
        & (x <= X_DAM_MAX + Z_TOL)
        & (y >= Y_DAM_MIN - Z_TOL)
        & (y <= Y_DAM_MAX + Z_TOL)
    )

    outside_dam_mask = ~dam_box_mask

    initial_positive_dam_mask = (
        dam_box_mask
        & (z_initial_geometry > Z_MIN + Z_TOL)
    )

    print(
        "Nodes in dam activation box :",
        int(np.sum(dam_box_mask))
    )

    print(
        "Initially elevated dam nodes:",
        int(np.sum(initial_positive_dam_mask))
    )

    print(
        "Nodes outside dam box        :",
        int(np.sum(outside_dam_mask))
    )

    print()


    # --------------------------------------------------------
    # GLOBAL STATUS FLAGS
    # --------------------------------------------------------

    passed = True

    previous_max_lowering = -1.0

    initial_result_bottom = None


    # --------------------------------------------------------
    # CHECK SELECTED TIMES
    # --------------------------------------------------------

    print_separator()

    for target_time in CHECK_TIMES:

        record = nearest_record(times, target_time)

        actual_time = times[record]

        z = np.asarray(
            res.get_data_value(bottom_res_name, record),
            dtype=float
        )

        if record == nearest_record(times, 0.0):
            initial_result_bottom = z.copy()

        lowering = z_initial_geometry - z

        raising = z - z_initial_geometry


        # ----------------------------------------------------
        # BASIC FINITE CHECK
        # ----------------------------------------------------

        finite_ok = np.all(np.isfinite(z))

        if not finite_ok:
            passed = False


        # ----------------------------------------------------
        # INITIAL / CURRENT BED STATISTICS
        # ----------------------------------------------------

        min_z = float(np.nanmin(z))
        max_z = float(np.nanmax(z))

        max_lowering = float(np.nanmax(lowering))
        max_raising = float(np.nanmax(raising))

        eroded_mask = lowering > CHANGE_TOL

        raised_mask = raising > CHANGE_TOL

        eroded_nodes = int(np.sum(eroded_mask))
        raised_nodes = int(np.sum(raised_mask))


        # ----------------------------------------------------
        # DAM-REGION CHANGE
        # ----------------------------------------------------

        dam_lowering = lowering[dam_box_mask]

        if dam_lowering.size > 0:
            max_dam_lowering = float(
                np.nanmax(dam_lowering)
            )
        else:
            max_dam_lowering = 0.0


        # ----------------------------------------------------
        # OUTSIDE-DAM CHANGE
        # ----------------------------------------------------

        outside_lowering = lowering[outside_dam_mask]
        outside_raising = raising[outside_dam_mask]

        if outside_lowering.size > 0:

            max_outside_lowering = float(
                np.nanmax(outside_lowering)
            )

            max_outside_raising = float(
                np.nanmax(outside_raising)
            )

            changed_outside = int(
                np.sum(
                    np.abs(
                        z[outside_dam_mask]
                        - z_initial_geometry[outside_dam_mask]
                    ) > CHANGE_TOL
                )
            )

        else:

            max_outside_lowering = 0.0
            max_outside_raising = 0.0
            changed_outside = 0


        # ----------------------------------------------------
        # BELOW FLUME FLOOR
        # ----------------------------------------------------

        below_floor_mask = z < Z_MIN - CHANGE_TOL

        below_floor_nodes = int(
            np.sum(below_floor_mask)
        )


        # ----------------------------------------------------
        # TEMPORAL ACCUMULATION
        # ----------------------------------------------------

        monotonic_max_erosion = (
            max_lowering + CHANGE_TOL
            >= previous_max_lowering
        )

        if not monotonic_max_erosion:
            passed = False

        previous_max_lowering = max(
            previous_max_lowering,
            max_lowering
        )


        # ----------------------------------------------------
        # STATUS CHECKS
        # ----------------------------------------------------

        if max_raising > CHANGE_TOL:
            passed = False

        if changed_outside > 0:
            passed = False

        if below_floor_nodes > 0:
            passed = False


        # ----------------------------------------------------
        # PRINT
        # ----------------------------------------------------

        print(
            f"Target time             : {target_time:8.2f} s"
        )

        print(
            f"Actual output time      : {actual_time:8.6f} s"
        )

        print(
            f"Finite bottom           : {finite_ok}"
        )

        print(
            f"Minimum bottom          : {min_z:.10f} m"
        )

        print(
            f"Maximum bottom          : {max_z:.10f} m"
        )

        print(
            f"Maximum lowering        : {max_lowering:.10f} m"
        )

        print(
            f"Maximum dam lowering    : {max_dam_lowering:.10f} m"
        )

        print(
            f"Eroded nodes            : {eroded_nodes}"
        )

        print(
            f"Maximum bed raising     : {max_raising:.10e} m"
        )

        print(
            f"Raised nodes            : {raised_nodes}"
        )

        print(
            f"Changed outside dam     : {changed_outside}"
        )

        print(
            f"Max outside lowering    : "
            f"{max_outside_lowering:.10e} m"
        )

        print(
            f"Max outside raising     : "
            f"{max_outside_raising:.10e} m"
        )

        print(
            f"Nodes below Z=0         : {below_floor_nodes}"
        )

        print(
            f"Max erosion nondecrease : "
            f"{monotonic_max_erosion}"
        )

        print_separator()


    # --------------------------------------------------------
    # INITIAL BOTTOM CONSISTENCY
    # --------------------------------------------------------

    if initial_result_bottom is None:
        raise RuntimeError(
            "Could not obtain initial result bottom."
        )

    initial_difference = (
        initial_result_bottom - z_initial_geometry
    )

    max_initial_abs_difference = float(
        np.nanmax(np.abs(initial_difference))
    )

    initial_consistency_ok = (
        max_initial_abs_difference <= CHANGE_TOL
    )

    if not initial_consistency_ok:
        passed = False


    # --------------------------------------------------------
    # FINAL BED
    # --------------------------------------------------------

    final_record = nearest_record(times, 30.0)

    z_final = np.asarray(
        res.get_data_value(
            bottom_res_name,
            final_record
        ),
        dtype=float
    )

    final_lowering = (
        z_initial_geometry - z_final
    )

    final_eroded_mask = (
        final_lowering > CHANGE_TOL
    )

    final_eroded_dam_mask = (
        final_eroded_mask & dam_box_mask
    )

    final_eroded_outside_mask = (
        final_eroded_mask & outside_dam_mask
    )


    # --------------------------------------------------------
    # FINAL SUMMARY
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL REFERENCE-RUN SUMMARY")
    print("=" * 78)

    print(
        "Initial bottom max |difference| : "
        f"{max_initial_abs_difference:.10e} m"
    )

    print(
        "Initial bottom consistency       : "
        f"{initial_consistency_ok}"
    )

    print(
        "Final maximum lowering           : "
        f"{np.max(final_lowering):.10f} m"
    )

    print(
        "Final eroded nodes               : "
        f"{np.sum(final_eroded_mask)}"
    )

    print(
        "Final eroded nodes in dam box    : "
        f"{np.sum(final_eroded_dam_mask)}"
    )

    print(
        "Final eroded nodes outside dam   : "
        f"{np.sum(final_eroded_outside_mask)}"
    )

    print(
        "Final minimum bottom             : "
        f"{np.min(z_final):.10f} m"
    )

    print(
        "Final maximum bottom             : "
        f"{np.max(z_final):.10f} m"
    )


    # --------------------------------------------------------
    # FINAL DECISION
    # --------------------------------------------------------

    if passed:

        print()
        print(
            "STATUS: B2 EROSION REFERENCE "
            "IMPLEMENTATION CHECK PASSED"
        )

        print()
        print(
            "NOTE: This confirms numerical / implementation "
            "behaviour only."
        )

        print(
            "      It does NOT mean that tau_c and k_d are "
            "calibrated."
        )

    else:

        print()
        print(
            "STATUS: B2 EROSION REFERENCE "
            "IMPLEMENTATION CHECK FAILED"
        )

        print(
            "Do not proceed to parameter calibration."
        )


finally:

    if geo is not None:
        geo.close()

    if res is not None:
        res.close()
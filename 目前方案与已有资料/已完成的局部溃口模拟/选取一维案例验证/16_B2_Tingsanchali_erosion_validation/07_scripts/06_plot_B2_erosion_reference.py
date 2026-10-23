import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# B2 Tingsanchali erosion reference-run visualization
#
# TELEMAC 9.1 compatible version
#
# Outputs:
#   1. B2_erosion_ref_longitudinal_profile.png
#   2. B2_fixedbed_vs_erosion_reservoir_level.png
#   3. B2_fixedbed_vs_erosion_crest_discharge.png
# ============================================================


# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

ROOT = (
    r"F:\Telemac_Projects\Tuanzho_New"
    r"\16_B2_Tingsanchali_erosion_validation"
)

FIXED_FILE = os.path.join(
    ROOT,
    "05_run",
    "B2_Tingsanchali_fixedbed_res.slf"
)

EROSION_FILE = os.path.join(
    ROOT,
    "05_run",
    "B2_Tingsanchali_erosion_ref_res.slf"
)

FIG_DIR = os.path.join(
    ROOT,
    "08_figures"
)

os.makedirs(FIG_DIR, exist_ok=True)


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

PROFILE_TIMES = [
    0.0,
    5.0,
    10.0,
    20.0,
    30.0,
]

TIME_TOL = 0.51

H_DRY = 1.0e-4

Y_CENTER = 0.50

PROFILE_START = [-17.50, Y_CENTER]
PROFILE_END = [17.50, Y_CENTER]

RES_X_MIN = -10.0
RES_X_MAX = -5.0

Q_SECTION_X = 0.0


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

def find_variable_name(varnames, candidates):

    for candidate in candidates:

        for name in varnames:

            clean = name.strip().upper()

            if clean == candidate.upper():
                return name

    return None


def nearest_record(times, target):

    times = np.asarray(times, dtype=float)

    idx = int(
        np.argmin(
            np.abs(times - target)
        )
    )

    if abs(times[idx] - target) > TIME_TOL:

        raise RuntimeError(
            f"No output record sufficiently close "
            f"to t={target:.2f} s. "
            f"Nearest = {times[idx]:.6f} s."
        )

    return idx


def get_polyline_timeseries(
    telemac_file,
    variable_name,
    polyline,
    poly_number
):

    """
    TELEMAC 9.1:
    get_timeseries_on_polyline() returns:

        points,
        curvilinear_abscissa,
        values

    values normally has shape:
        (number_of_profile_points, number_of_times)
    """

    result = telemac_file.get_timeseries_on_polyline(
        variable_name,
        polyline,
        poly_number
    )

    if len(result) != 3:

        raise RuntimeError(
            "Unexpected TELEMAC polyline return format. "
            f"Number of returned objects = {len(result)}"
        )

    points = np.asarray(
        result[0],
        dtype=float
    )

    curvilinear_abscissa = np.asarray(
        result[1],
        dtype=float
    )

    values = np.asarray(
        result[2],
        dtype=float
    )

    return (
        points,
        curvilinear_abscissa,
        values
    )


def extract_record_from_polyline_values(
    values,
    n_times,
    record
):

    """
    Handle TELEMAC polyline array orientation defensively.
    """

    values = np.asarray(
        values,
        dtype=float
    )

    if values.ndim != 2:

        raise RuntimeError(
            "Unexpected polyline values dimension: "
            f"{values.shape}"
        )

    if values.shape[1] == n_times:

        return values[:, record]

    if values.shape[0] == n_times:

        return values[record, :]

    raise RuntimeError(
        "Unexpected polyline values shape: "
        f"{values.shape}; n_times={n_times}"
    )


def get_centerline_profile(
    telemac_file,
    variable_name,
    record
):

    """
    Extract longitudinal profile along y = 0.5 m.
    """

    polyline = [
        PROFILE_START,
        PROFILE_END,
    ]

    # 35 m / 0.05 m = 700 intervals
    poly_number = [701]

    points, _, values = get_polyline_timeseries(
        telemac_file,
        variable_name,
        polyline,
        poly_number
    )

    profile = extract_record_from_polyline_values(
        values,
        len(telemac_file.times),
        record
    )

    x_profile = points[:, 0]

    return (
        x_profile,
        np.asarray(profile, dtype=float)
    )


def reservoir_level_series(
    telemac_file,
    free_surface_name
):

    """
    Mean free-surface elevation over:
        -10 <= x <= -5 m
    """

    x = np.asarray(
        telemac_file.meshx,
        dtype=float
    )

    mask = (
        (x >= RES_X_MIN)
        & (x <= RES_X_MAX)
    )

    if np.sum(mask) == 0:

        raise RuntimeError(
            "No nodes found in reservoir sampling region."
        )

    eta_series = np.zeros(
        len(telemac_file.times),
        dtype=float
    )

    for record in range(
        len(telemac_file.times)
    ):

        eta = np.asarray(
            telemac_file.get_data_value(
                free_surface_name,
                record
            ),
            dtype=float
        )

        eta_series[record] = float(
            np.mean(
                eta[mask]
            )
        )

    return eta_series


def discharge_series(
    telemac_file,
    h_name,
    u_name
):

    """
    Compute:
        Q = integral(H * U dy)

    across x = 0 from y = 0 to 1 m.
    """

    polyline = [
        [Q_SECTION_X, 0.0],
        [Q_SECTION_X, 1.0],
    ]

    # Use dense interpolation across the 1-m-wide flume.
    poly_number = [101]

    points_h, s_h, h_values = get_polyline_timeseries(
        telemac_file,
        h_name,
        polyline,
        poly_number
    )

    points_u, s_u, u_values = get_polyline_timeseries(
        telemac_file,
        u_name,
        polyline,
        poly_number
    )

    if len(points_h) != len(points_u):

        raise RuntimeError(
            "H and U section point counts differ."
        )

    y_section = points_h[:, 1]

    order = np.argsort(
        y_section
    )

    y_section = y_section[order]

    n_times = len(
        telemac_file.times
    )

    q = np.zeros(
        n_times,
        dtype=float
    )

    for record in range(
        n_times
    ):

        h = extract_record_from_polyline_values(
            h_values,
            n_times,
            record
        )

        u = extract_record_from_polyline_values(
            u_values,
            n_times,
            record
        )

        h = h[order]
        u = u[order]

        q_density = h * u

        q[record] = np.trapz(
            q_density,
            y_section
        )

    return q


# ------------------------------------------------------------
# FILE CHECK
# ------------------------------------------------------------

print()
print("=" * 78)
print("B2 EROSION REFERENCE VISUALIZATION")
print("=" * 78)

for path in [
    FIXED_FILE,
    EROSION_FILE,
]:

    if not os.path.isfile(path):

        print("ERROR: File not found:")
        print(path)

        sys.exit(1)

print("Fixed-bed result :", FIXED_FILE)
print("Erosion result   :", EROSION_FILE)
print("Figure directory :", FIG_DIR)
print()


# ------------------------------------------------------------
# OPEN FILES
# ------------------------------------------------------------

fixed = None
erosion = None

try:

    fixed = TelemacFile(
        FIXED_FILE
    )

    erosion = TelemacFile(
        EROSION_FILE
    )


    # --------------------------------------------------------
    # VARIABLES
    # --------------------------------------------------------

    fixed_h = find_variable_name(
        fixed.varnames,
        ["WATER DEPTH", "H"]
    )

    fixed_u = find_variable_name(
        fixed.varnames,
        ["VELOCITY U", "U"]
    )

    fixed_s = find_variable_name(
        fixed.varnames,
        ["FREE SURFACE", "S"]
    )

    fixed_b = find_variable_name(
        fixed.varnames,
        ["BOTTOM", "B"]
    )


    erosion_h = find_variable_name(
        erosion.varnames,
        ["WATER DEPTH", "H"]
    )

    erosion_u = find_variable_name(
        erosion.varnames,
        ["VELOCITY U", "U"]
    )

    erosion_s = find_variable_name(
        erosion.varnames,
        ["FREE SURFACE", "S"]
    )

    erosion_b = find_variable_name(
        erosion.varnames,
        ["BOTTOM", "B"]
    )


    required = {
        "fixed WATER DEPTH": fixed_h,
        "fixed VELOCITY U": fixed_u,
        "fixed FREE SURFACE": fixed_s,
        "fixed BOTTOM": fixed_b,
        "erosion WATER DEPTH": erosion_h,
        "erosion VELOCITY U": erosion_u,
        "erosion FREE SURFACE": erosion_s,
        "erosion BOTTOM": erosion_b,
    }

    for label, value in required.items():

        if value is None:

            raise RuntimeError(
                f"Required variable not found: {label}"
            )


    print("Variables detected successfully.")
    print()


    # ========================================================
    # FIGURE 1
    #
    # LONGITUDINAL BED + WATER-SURFACE EVOLUTION
    # ========================================================

    print(
        "Creating longitudinal erosion-profile figure..."
    )

    fig, ax = plt.subplots(
        figsize=(12, 7)
    )


    # --------------------------------------------------------
    # Initial bed
    # --------------------------------------------------------

    record0 = nearest_record(
        erosion.times,
        0.0
    )

    x_profile, bed0 = get_centerline_profile(
        erosion,
        erosion_b,
        record0
    )

    ax.plot(
        x_profile,
        bed0,
        linewidth=2.5,
        label="Initial bed"
    )


    # --------------------------------------------------------
    # Bed + water surface at selected times
    # --------------------------------------------------------

    for target_time in PROFILE_TIMES:

        record = nearest_record(
            erosion.times,
            target_time
        )

        actual_time = float(
            erosion.times[record]
        )


        x_bed, bed = get_centerline_profile(
            erosion,
            erosion_b,
            record
        )

        x_eta, eta = get_centerline_profile(
            erosion,
            erosion_s,
            record
        )

        x_h, depth = get_centerline_profile(
            erosion,
            erosion_h,
            record
        )


        eta_plot = eta.copy()

        eta_plot[
            depth <= H_DRY
        ] = np.nan


        if target_time > 0.0:

            ax.plot(
                x_bed,
                bed,
                linewidth=1.7,
                label=(
                    f"Bed t={actual_time:.0f} s"
                )
            )


        ax.plot(
            x_eta,
            eta_plot,
            linestyle="--",
            linewidth=1.3,
            label=(
                f"Water surface "
                f"t={actual_time:.0f} s"
            )
        )


    # --------------------------------------------------------
    # Important dam locations
    # --------------------------------------------------------

    key_positions = [
        -2.40,
        0.00,
        0.30,
        4.30,
    ]

    for xpos in key_positions:

        ax.axvline(
            xpos,
            linestyle=":",
            linewidth=0.9
        )


    # Focus on the dam region.
    # Full flume is 35 m long, but the erosion morphology
    # is much clearer at this scale.

    ax.set_xlim(
        -3.5,
        5.5
    )

    ax.set_xlabel(
        "Longitudinal distance x (m)"
    )

    ax.set_ylabel(
        "Elevation (m)"
    )

    ax.set_title(
        "B-2 erosion reference run: "
        "longitudinal bed and water-surface evolution"
    )

    ax.grid(
        True,
        alpha=0.25
    )

    ax.legend(
        fontsize=8,
        ncol=2
    )

    fig.tight_layout()


    output_1 = os.path.join(
        FIG_DIR,
        "B2_erosion_ref_longitudinal_profile.png"
    )

    fig.savefig(
        output_1,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print("Saved:")
    print(output_1)
    print()


    # ========================================================
    # FIGURE 2
    #
    # RESERVOIR LEVEL
    # FIXED BED vs ERODIBLE BED
    # ========================================================

    print(
        "Creating reservoir-level comparison..."
    )

    fixed_eta = reservoir_level_series(
        fixed,
        fixed_s
    )

    erosion_eta = reservoir_level_series(
        erosion,
        erosion_s
    )


    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.plot(
        fixed.times,
        fixed_eta,
        linewidth=2.0,
        label="Fixed bed"
    )

    ax.plot(
        erosion.times,
        erosion_eta,
        linewidth=2.0,
        label="Erodible bed"
    )

    ax.axhline(
        0.85,
        linestyle=":",
        linewidth=1.0,
        label="Initial reservoir level"
    )

    ax.set_xlabel(
        "Time (s)"
    )

    ax.set_ylabel(
        "Mean reservoir free-surface elevation (m)"
    )

    ax.set_title(
        "B-2 reservoir level: "
        "fixed bed vs erodible bed"
    )

    ax.grid(
        True,
        alpha=0.25
    )

    ax.legend()

    fig.tight_layout()


    output_2 = os.path.join(
        FIG_DIR,
        "B2_fixedbed_vs_erosion_reservoir_level.png"
    )

    fig.savefig(
        output_2,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print("Saved:")
    print(output_2)
    print()


    # ========================================================
    # FIGURE 3
    #
    # CREST-FRONT DISCHARGE
    # FIXED BED vs ERODIBLE BED
    # ========================================================

    print(
        "Computing crest-front discharge series..."
    )

    fixed_q = discharge_series(
        fixed,
        fixed_h,
        fixed_u
    )

    erosion_q = discharge_series(
        erosion,
        erosion_h,
        erosion_u
    )


    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.plot(
        fixed.times,
        fixed_q,
        linewidth=2.0,
        label="Fixed bed"
    )

    ax.plot(
        erosion.times,
        erosion_q,
        linewidth=2.0,
        label="Erodible bed"
    )

    ax.set_xlabel(
        "Time (s)"
    )

    ax.set_ylabel(
        "Discharge Q at x = 0 (m$^3$/s)"
    )

    ax.set_title(
        "B-2 crest-front discharge: "
        "fixed bed vs erodible bed"
    )

    ax.grid(
        True,
        alpha=0.25
    )

    ax.legend()

    fig.tight_layout()


    output_3 = os.path.join(
        FIG_DIR,
        "B2_fixedbed_vs_erosion_crest_discharge.png"
    )

    fig.savefig(
        output_3,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print("Saved:")
    print(output_3)
    print()


    # ========================================================
    # NUMERICAL SUMMARY
    # ========================================================

    print("=" * 78)
    print("VISUALIZATION SUMMARY")
    print("=" * 78)

    print(
        "Fixed-bed reservoir level:"
    )

    print(
        f"  t=0 s  : "
        f"{fixed_eta[0]:.6f} m"
    )

    print(
        f"  t=30 s : "
        f"{fixed_eta[-1]:.6f} m"
    )

    print()

    print(
        "Erosion-run reservoir level:"
    )

    print(
        f"  t=0 s  : "
        f"{erosion_eta[0]:.6f} m"
    )

    print(
        f"  t=30 s : "
        f"{erosion_eta[-1]:.6f} m"
    )

    print()

    print(
        "Difference in reservoir level "
        "at t=30 s:"
    )

    print(
        f"  erosion - fixed = "
        f"{erosion_eta[-1] - fixed_eta[-1]:.8f} m"
    )

    print()

    print(
        "Fixed-bed crest discharge:"
    )

    print(
        f"  max    : "
        f"{np.max(fixed_q):.8f} m3/s"
    )

    print(
        f"  t=30 s : "
        f"{fixed_q[-1]:.8f} m3/s"
    )

    print()

    print(
        "Erosion-run crest discharge:"
    )

    print(
        f"  max    : "
        f"{np.max(erosion_q):.8f} m3/s"
    )

    print(
        f"  t=30 s : "
        f"{erosion_q[-1]:.8f} m3/s"
    )

    print()

    if np.max(fixed_q) != 0.0:

        qmax_change = (
            np.max(erosion_q)
            - np.max(fixed_q)
        ) / np.max(fixed_q) * 100.0

        print(
            "Peak-discharge change:"
        )

        print(
            f"  {qmax_change:.3f} %"
        )

        print()


    print(
        "Created figures:"
    )

    print(
        "  1.",
        os.path.basename(output_1)
    )

    print(
        "  2.",
        os.path.basename(output_2)
    )

    print(
        "  3.",
        os.path.basename(output_3)
    )

    print()

    print(
        "STATUS: B2 EROSION REFERENCE FIGURES CREATED"
    )

    print("=" * 78)


finally:

    if fixed is not None:
        fixed.close()

    if erosion is not None:
        erosion.close()
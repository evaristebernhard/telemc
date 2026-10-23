# -*- coding: utf-8 -*-

"""
TUANZHOU L2-B1 600 s erosion-evolution visualisation

Purpose
-------
Visualise cumulative vertical erosion from the completed
L2-B1 600 s TELEMAC result.

The plotted quantity is

    erosion_depth(x,y,t) = BOTTOM(x,y,0) - BOTTOM(x,y,t)

Therefore:

    erosion_depth = 0     -> no erosion
    erosion_depth > 0     -> cumulative erosion

This script creates:

1. Four individual cumulative-erosion maps:
       t = 60, 120, 300, 600 s

2. One detailed cumulative-erosion map at 600 s.

3. One four-panel comparison figure:
       60 | 120 | 300 | 600 s

4. One B1 time-series figure of maximum cumulative erosion depth.

5. One B0/B1 comparison figure of maximum cumulative erosion depth.

No TELEMAC model files are modified.
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


# ======================================================================
# PATHS
# ======================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

RUN_DIR = PROJECT_DIR / "05_run"

OUTPUT_DIR = (
    PROJECT_DIR
    / "06_results"
    / "B1_600s_erosion_figures"
)

B1_RESULT_FILE = (
    RUN_DIR
    / "local_breach_L2_B1_600s_res.slf"
)

B0_RESULT_FILE = Path(
    r"F:\Telemac_Projects\Tuanzho_New"
    r"\14_local_breach_L2_B0_vertical_erosion"
    r"\09_B0_600s_test"
    r"\05_run"
    r"\local_breach_L2_B0_600s_res.slf"
)


# ======================================================================
# VISUALISATION PARAMETERS
# ======================================================================

MAP_TIMES = [
    60.0,
    120.0,
    300.0,
    600.0,
]

X_MIN_PLOT = -20.0
X_MAX_PLOT = 20.0

Y_MIN_PLOT = -20.0
Y_MAX_PLOT = 25.0

INITIAL_BREACH_LEFT = -10.0
INITIAL_BREACH_RIGHT = 10.0

EROSION_PLOT_TOL = 1.0e-6

N_LEVELS = 31


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
    print(
        "Run this script inside the configured "
        "openTELEMAC environment."
    )
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


def get_triangles(res):

    candidates = [
        "ikle2",
        "ikle",
    ]

    for name in candidates:

        if hasattr(res, name):

            triangles = np.asarray(
                getattr(res, name)
            )

            if triangles.ndim == 2:

                if triangles.shape[1] >= 3:

                    triangles = triangles[:, :3]

                    return triangles.astype(int)

    print()
    print("ERROR:")
    print(
        "Could not obtain triangular connectivity "
        "from the TELEMAC result."
    )
    print()

    sys.exit(1)


def cumulative_erosion(bottom0, bottom):

    erosion = bottom0 - bottom

    erosion[
        np.abs(erosion) < EROSION_PLOT_TOL
    ] = 0.0

    return erosion


def make_zoom_mask(x, y):

    return (
        (x >= X_MIN_PLOT)
        & (x <= X_MAX_PLOT)
        & (y >= Y_MIN_PLOT)
        & (y <= Y_MAX_PLOT)
    )


def make_zoom_triangulation(
    x,
    y,
    triangles,
):

    triangulation = mtri.Triangulation(
        x,
        y,
        triangles,
    )

    tri_x = x[triangles]
    tri_y = y[triangles]

    tri_centroid_x = np.mean(
        tri_x,
        axis=1,
    )

    tri_centroid_y = np.mean(
        tri_y,
        axis=1,
    )

    outside = (
        (tri_centroid_x < X_MIN_PLOT)
        | (tri_centroid_x > X_MAX_PLOT)
        | (tri_centroid_y < Y_MIN_PLOT)
        | (tri_centroid_y > Y_MAX_PLOT)
    )

    triangulation.set_mask(
        outside
    )

    return triangulation


def add_reference_lines(ax):

    ax.axvline(
        INITIAL_BREACH_LEFT,
        linestyle="--",
        linewidth=1.0,
    )

    ax.axvline(
        INITIAL_BREACH_RIGHT,
        linestyle="--",
        linewidth=1.0,
    )

    ax.axhline(
        0.0,
        linestyle=":",
        linewidth=1.0,
    )


def format_spatial_axis(ax):

    ax.set_xlim(
        X_MIN_PLOT,
        X_MAX_PLOT,
    )

    ax.set_ylim(
        Y_MIN_PLOT,
        Y_MAX_PLOT,
    )

    ax.set_aspect(
        "equal",
        adjustable="box",
    )

    ax.set_xlabel(
        "x (m)"
    )

    ax.set_ylabel(
        "y (m)"
    )

    add_reference_lines(
        ax
    )


def strongest_point(
    erosion,
    zoom_mask,
):

    local_indices = np.where(
        zoom_mask
    )[0]

    return int(
        local_indices[
            np.argmax(
                erosion[local_indices]
            )
        ]
    )


def calculate_max_series(
    res,
    bottom0,
    zoom_mask,
):

    times = np.asarray(
        res.times,
        dtype=float,
    )

    series = []

    for rec in range(len(times)):

        bottom = get_bottom(
            res,
            rec,
        )

        erosion = cumulative_erosion(
            bottom0,
            bottom,
        )

        series.append(
            float(
                np.max(
                    erosion[zoom_mask]
                )
            )
        )

    return np.asarray(
        series,
        dtype=float,
    )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print()
    print("=" * 78)
    print(
        "TUANZHOU L2-B1 600 s EROSION EVOLUTION VISUALISATION"
    )
    print("=" * 78)
    print()

    print("B1 result file:")
    print(B1_RESULT_FILE)
    print()

    if not B1_RESULT_FILE.exists():

        print("ERROR: B1 result file does not exist.")
        sys.exit(1)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    res = TelemacFile(
        str(B1_RESULT_FILE)
    )

    try:

        # ==============================================================
        # READ BASIC DATA
        # ==============================================================

        x = np.asarray(
            res.meshx,
            dtype=float,
        )

        y = np.asarray(
            res.meshy,
            dtype=float,
        )

        times = np.asarray(
            res.times,
            dtype=float,
        )

        triangles = get_triangles(
            res
        )

        print(
            f"Number of nodes      : {len(x)}"
        )

        print(
            f"Number of triangles  : {len(triangles)}"
        )

        print(
            f"Number of records    : {len(times)}"
        )

        print(
            f"Result time range    : "
            f"{times[0]:.3f} - "
            f"{times[-1]:.3f} s"
        )

        print()

        # ==============================================================
        # INITIAL BOTTOM
        # ==============================================================

        bottom0 = get_bottom(
            res,
            0,
        )

        # ==============================================================
        # TRIANGULATION AND MASK
        # ==============================================================

        triangulation = make_zoom_triangulation(
            x,
            y,
            triangles,
        )

        zoom_mask = make_zoom_mask(
            x,
            y,
        )

        # ==============================================================
        # READ COMPARISON TIMES
        # ==============================================================

        map_data = {}

        common_max_erosion = 0.0

        for target_time in MAP_TIMES:

            rec = nearest_record(
                times,
                target_time,
            )

            actual_time = float(
                times[rec]
            )

            bottom = get_bottom(
                res,
                rec,
            )

            erosion = cumulative_erosion(
                bottom0,
                bottom,
            )

            max_zoom_erosion = float(
                np.max(
                    erosion[zoom_mask]
                )
            )

            common_max_erosion = max(
                common_max_erosion,
                max_zoom_erosion,
            )

            map_data[target_time] = {
                "record": rec,
                "actual_time": actual_time,
                "erosion": erosion,
                "max_erosion": max_zoom_erosion,
            }

        if common_max_erosion <= 0.0:

            print("ERROR:")
            print(
                "No positive cumulative erosion "
                "was detected."
            )

            sys.exit(1)

        levels = np.linspace(
            0.0,
            common_max_erosion,
            N_LEVELS,
        )

        print("Common map colour scale:")
        print(
            f"0.000000 - "
            f"{common_max_erosion:.9f} m"
        )
        print()

        # ==============================================================
        # FIGURES 1-4
        # ==============================================================

        for target_time in MAP_TIMES:

            data = map_data[
                target_time
            ]

            actual_time = data[
                "actual_time"
            ]

            erosion = data[
                "erosion"
            ]

            strongest = strongest_point(
                erosion,
                zoom_mask,
            )

            fig, ax = plt.subplots(
                figsize=(8.0, 7.0)
            )

            contour = ax.tricontourf(
                triangulation,
                erosion,
                levels=levels,
                extend="max",
            )

            cbar = fig.colorbar(
                contour,
                ax=ax,
                pad=0.02,
            )

            cbar.set_label(
                "Cumulative erosion depth (m)"
            )

            format_spatial_axis(
                ax
            )

            ax.set_title(
                "L2-B1 cumulative vertical erosion\n"
                f"t = {actual_time:.1f} s"
            )

            ax.plot(
                x[strongest],
                y[strongest],
                marker="x",
                markersize=9,
                markeredgewidth=2,
            )

            ax.text(
                x[strongest] + 0.7,
                y[strongest] + 0.7,
                (
                    f"max = "
                    f"{erosion[strongest]:.4f} m"
                ),
                fontsize=9,
            )

            fig.tight_layout()

            output_file = (
                OUTPUT_DIR
                / (
                    "B1_600s_cumulative_erosion_"
                    f"{int(round(actual_time)):03d}s.png"
                )
            )

            fig.savefig(
                output_file,
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(
                fig
            )

            print(
                f"Saved: {output_file.name}"
            )

        # ==============================================================
        # FIGURE 5
        # Detailed 600 s erosion map
        # ==============================================================

        final_target_time = MAP_TIMES[-1]

        final_data = map_data[
            final_target_time
        ]

        final_time = final_data[
            "actual_time"
        ]

        final_erosion = final_data[
            "erosion"
        ]

        final_max = float(
            np.max(
                final_erosion[
                    zoom_mask
                ]
            )
        )

        detailed_levels = np.linspace(
            0.0,
            final_max,
            41,
        )

        fig, ax = plt.subplots(
            figsize=(8.5, 7.5)
        )

        contour = ax.tricontourf(
            triangulation,
            final_erosion,
            levels=detailed_levels,
            extend="max",
        )

        contour_values = [
            0.01,
            0.05,
            0.10,
            0.15,
            0.20,
            0.25,
        ]

        contour_values = [
            value
            for value in contour_values
            if value < final_max
        ]

        if len(contour_values) > 0:

            line_contours = ax.tricontour(
                triangulation,
                final_erosion,
                levels=contour_values,
                linewidths=0.8,
            )

            ax.clabel(
                line_contours,
                inline=True,
                fontsize=8,
                fmt="%.2f m",
            )

        cbar = fig.colorbar(
            contour,
            ax=ax,
            pad=0.02,
        )

        cbar.set_label(
            "Cumulative erosion depth (m)"
        )

        format_spatial_axis(
            ax
        )

        ax.set_title(
            "L2-B1 detailed cumulative erosion\n"
            f"t = {final_time:.1f} s"
        )

        strongest_final = strongest_point(
            final_erosion,
            zoom_mask,
        )

        ax.plot(
            x[strongest_final],
            y[strongest_final],
            marker="x",
            markersize=10,
            markeredgewidth=2,
        )

        ax.text(
            x[strongest_final] + 0.7,
            y[strongest_final] + 0.7,
            (
                "Maximum erosion\n"
                f"{final_erosion[strongest_final]:.5f} m"
            ),
            fontsize=9,
        )

        fig.tight_layout()

        detailed_output = (
            OUTPUT_DIR
            / "B1_600s_detailed_erosion_600s.png"
        )

        fig.savefig(
            detailed_output,
            dpi=300,
            bbox_inches="tight",
        )

        plt.close(
            fig
        )

        print(
            f"Saved: {detailed_output.name}"
        )

        # ==============================================================
        # FIGURE 6
        # Four-panel comparison
        # ==============================================================

        fig, axes = plt.subplots(
            1,
            4,
            figsize=(19.0, 6.0),
            sharex=True,
            sharey=True,
        )

        comparison_contour = None

        for ax, target_time in zip(
            axes,
            MAP_TIMES,
        ):

            data = map_data[
                target_time
            ]

            actual_time = data[
                "actual_time"
            ]

            erosion = data[
                "erosion"
            ]

            comparison_contour = ax.tricontourf(
                triangulation,
                erosion,
                levels=levels,
                extend="max",
            )

            format_spatial_axis(
                ax
            )

            ax.set_title(
                f"t = {actual_time:.0f} s"
            )

            strongest = strongest_point(
                erosion,
                zoom_mask,
            )

            ax.plot(
                x[strongest],
                y[strongest],
                marker="x",
                markersize=7,
                markeredgewidth=1.5,
            )

            ax.text(
                0.03,
                0.97,
                (
                    "max = "
                    f"{erosion[strongest]:.3f} m"
                ),
                transform=ax.transAxes,
                verticalalignment="top",
                fontsize=9,
            )

        for ax in axes[1:]:
            ax.set_ylabel("")

        cbar = fig.colorbar(
            comparison_contour,
            ax=axes,
            pad=0.02,
            fraction=0.025,
        )

        cbar.set_label(
            "Cumulative erosion depth (m)"
        )

        fig.suptitle(
            "L2-B1 cumulative vertical erosion evolution",
            fontsize=14,
        )

        fig.subplots_adjust(
            left=0.055,
            right=0.91,
            bottom=0.12,
            top=0.86,
            wspace=0.18,
        )

        comparison_output = (
            OUTPUT_DIR
            / "B1_600s_erosion_comparison_060_120_300_600s.png"
        )

        fig.savefig(
            comparison_output,
            dpi=300,
            bbox_inches="tight",
        )

        plt.close(
            fig
        )

        print(
            f"Saved: {comparison_output.name}"
        )

        # ==============================================================
        # FIGURE 7
        # B1 maximum cumulative erosion versus time
        # ==============================================================

        max_erosion_series = calculate_max_series(
            res,
            bottom0,
            zoom_mask,
        )

        fig, ax = plt.subplots(
            figsize=(9.0, 5.8)
        )

        ax.plot(
            times,
            max_erosion_series,
            linewidth=1.8,
        )

        ax.set_xlabel(
            "Time (s)"
        )

        ax.set_ylabel(
            "Maximum cumulative erosion depth (m)"
        )

        ax.set_title(
            "L2-B1 maximum cumulative erosion evolution"
        )

        ax.grid(
            True,
            linewidth=0.5,
            alpha=0.35,
        )

        for target_time in MAP_TIMES:

            rec = nearest_record(
                times,
                target_time,
            )

            ax.plot(
                times[rec],
                max_erosion_series[rec],
                marker="o",
                markersize=5,
            )

            ax.annotate(
                (
                    f"{times[rec]:.0f} s\n"
                    f"{max_erosion_series[rec]:.3f} m"
                ),
                xy=(
                    times[rec],
                    max_erosion_series[rec],
                ),
                xytext=(6, 8),
                textcoords="offset points",
                fontsize=8,
            )

        fig.tight_layout()

        time_output = (
            OUTPUT_DIR
            / "B1_600s_maximum_erosion_vs_time.png"
        )

        fig.savefig(
            time_output,
            dpi=300,
            bbox_inches="tight",
        )

        plt.close(
            fig
        )

        print(
            f"Saved: {time_output.name}"
        )

        # ==============================================================
        # FIGURE 8
        # B0 versus B1 maximum erosion evolution
        # ==============================================================

        if B0_RESULT_FILE.exists():

            b0 = TelemacFile(
                str(B0_RESULT_FILE)
            )

            try:

                b0_x = np.asarray(
                    b0.meshx,
                    dtype=float,
                )

                b0_y = np.asarray(
                    b0.meshy,
                    dtype=float,
                )

                b0_times = np.asarray(
                    b0.times,
                    dtype=float,
                )

                if (
                    len(b0_x) != len(x)
                    or np.max(np.abs(b0_x - x)) > 1.0e-6
                    or np.max(np.abs(b0_y - y)) > 1.0e-6
                ):

                    print()
                    print(
                        "WARNING: B0/B1 meshes differ. "
                        "Skipping B0/B1 time-series comparison."
                    )

                else:

                    b0_bottom0 = get_bottom(
                        b0,
                        0,
                    )

                    b0_zoom_mask = make_zoom_mask(
                        b0_x,
                        b0_y,
                    )

                    b0_max_series = calculate_max_series(
                        b0,
                        b0_bottom0,
                        b0_zoom_mask,
                    )

                    fig, ax = plt.subplots(
                        figsize=(9.0, 5.8)
                    )

                    ax.plot(
                        b0_times,
                        b0_max_series,
                        linewidth=1.8,
                        label="L2-B0",
                    )

                    ax.plot(
                        times,
                        max_erosion_series,
                        linewidth=1.8,
                        label="L2-B1",
                    )

                    ax.set_xlabel(
                        "Time (s)"
                    )

                    ax.set_ylabel(
                        "Maximum cumulative erosion depth (m)"
                    )

                    ax.set_title(
                        "L2-B0 vs L2-B1 maximum cumulative erosion"
                    )

                    ax.grid(
                        True,
                        linewidth=0.5,
                        alpha=0.35,
                    )

                    ax.legend()

                    fig.tight_layout()

                    b0_b1_output = (
                        OUTPUT_DIR
                        / "B0_vs_B1_maximum_erosion_vs_time.png"
                    )

                    fig.savefig(
                        b0_b1_output,
                        dpi=300,
                        bbox_inches="tight",
                    )

                    plt.close(
                        fig
                    )

                    print(
                        f"Saved: {b0_b1_output.name}"
                    )

            finally:

                b0.close()

        else:

            print()
            print(
                "WARNING: B0 result file not found."
            )

            print(
                "B1 figures were created, but the "
                "B0/B1 comparison curve was skipped."
            )

        # ==============================================================
        # FINAL CONSOLE SUMMARY
        # ==============================================================

        print()
        print("=" * 78)
        print(
            "L2-B1 600 s EROSION VISUALISATION SUMMARY"
        )
        print("=" * 78)

        for target_time in MAP_TIMES:

            data = map_data[
                target_time
            ]

            print(
                f"t = "
                f"{data['actual_time']:6.1f} s"
                f" | max cumulative erosion = "
                f"{data['max_erosion']:.9e} m"
            )

        print("-" * 78)

        print(
            "Strongest erosion point at 600 s:"
        )

        print(
            f"x = {x[strongest_final]:.6f} m"
        )

        print(
            f"y = {y[strongest_final]:.6f} m"
        )

        print(
            f"erosion = "
            f"{final_erosion[strongest_final]:.9e} m"
        )

        print("-" * 78)

        print(
            "Figures saved to:"
        )

        print(
            OUTPUT_DIR
        )

        print("=" * 78)
        print()

    finally:

        res.close()


if __name__ == "__main__":
    main()
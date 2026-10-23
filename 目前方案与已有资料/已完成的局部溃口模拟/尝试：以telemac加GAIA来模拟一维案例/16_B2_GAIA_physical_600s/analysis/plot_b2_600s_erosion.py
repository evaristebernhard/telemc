from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize

from data_manip.extraction.telemac_file import TelemacFile

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "analysis"

GAIA = ROOT / "b2_physical_gaia_600s.slf"
T2D = ROOT / "b2_physical_t2d_600s.slf"

TIMES = [0, 60, 120, 300, 600]

DAM_X0 = 15.10
DAM_X1 = 21.80
CREST_X0 = 17.50
CREST_X1 = 17.80

plt.rcParams.update({
    "font.size": 10,
    "axes.grid": False,
    "savefig.dpi": 180
})

def find_variable(names, exact, contains=None):
    for wanted in exact:
        for name in names:
            if name.strip().upper() == wanted.upper():
                return name
    if contains:
        for wanted in contains:
            for name in names:
                if wanted.upper() in name.strip().upper():
                    return name
    return None

def nearest_index(times, target):
    return int(np.argmin(np.abs(times - target)))

def save_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

def main():
    if not GAIA.exists():
        raise FileNotFoundError(str(GAIA))

    OUT.mkdir(exist_ok=True)

    slf = TelemacFile(str(GAIA))
    try:
        names = [str(v).strip() for v in slf.varnames]

        print("GAIA variables:")
        for name in names:
            print("  ", name)

        bed_name = find_variable(
            names,
            ["BOTTOM", "BED ELEVATION", "BED LEVEL"],
            ["BOTTOM"]
        )

        evol_name = find_variable(
            names,
            ["EVOLUTION", "BED EVOLUTION"],
            ["EVOLUTION"]
        )

        if bed_name is None and evol_name is None:
            raise RuntimeError(
                "Neither BOTTOM nor EVOLUTION is available in GAIA SLF."
            )

        print("Bed variable:", bed_name)
        print("Evolution variable:", evol_name)

        x = np.asarray(slf.meshx, dtype=float)
        y = np.asarray(slf.meshy, dtype=float)
        triangles = np.asarray(slf.ikle2, dtype=int)
        times = np.asarray(slf.times, dtype=float)

        triang = mtri.Triangulation(x, y, triangles)

        if bed_name is not None:
            bed0 = np.asarray(
                slf.get_data_value(bed_name, 0), dtype=float
            )
        else:
            bed0 = None

        def get_change(it):
            if evol_name is not None:
                return np.asarray(
                    slf.get_data_value(evol_name, it), dtype=float
                )
            bed = np.asarray(
                slf.get_data_value(bed_name, it), dtype=float
            )
            return bed - bed0

        if bed0 is None:
            if T2D.exists():
                t2d = TelemacFile(str(T2D))
                try:
                    t2d_names = [
                        str(v).strip() for v in t2d.varnames
                    ]
                    t2d_bed = find_variable(
                        t2d_names,
                        ["BOTTOM", "BED ELEVATION", "BED LEVEL"],
                        ["BOTTOM"]
                    )
                    if t2d_bed is not None:
                        bed0 = np.asarray(
                            t2d.get_data_value(t2d_bed, 0),
                            dtype=float
                        )
                finally:
                    t2d.close()

        if bed0 is None:
            raise RuntimeError(
                "Initial bed elevation is unavailable."
            )

        dam_mask = (x >= DAM_X0) & (x <= DAM_X1)
        crest_mask = (x >= CREST_X0) & (x <= CREST_X1)

        if not np.any(crest_mask):
            raise RuntimeError("No nodes found on crest.")

        # Extract snapshots
        snapshots = {}
        for target in TIMES:
            it = nearest_index(times, target)
            snapshots[target] = {
                "actual_time": float(times[it]),
                "change": get_change(it)
            }

        # Figure 1: bed profiles
        # Select nodes close to flume centreline.
        yc = 0.5 * (float(y.min()) + float(y.max()))
        y_unique = np.unique(np.round(y, 6))
        dy = np.min(np.diff(y_unique)) if len(y_unique) > 1 else 0.05
        tol = max(float(dy) * 0.55, 0.025)

        center_ids = np.where(
            (np.abs(y - yc) <= tol) & dam_mask
        )[0]

        if len(center_ids) < 10:
            center_ids = np.where(dam_mask)[0]

        center_ids = center_ids[
            np.argsort(x[center_ids])
        ]

        fig, ax = plt.subplots(figsize=(11, 5.5))

        profile_rows = []

        for target in TIMES:
            dz = snapshots[target]["change"]
            z = bed0 + dz
            ax.plot(
                x[center_ids],
                z[center_ids],
                linewidth=1.5,
                label=f"{target} s"
            )
            for idx in center_ids:
                profile_rows.append([
                    target, x[idx], y[idx], z[idx], dz[idx]
                ])

        ax.axvspan(
            CREST_X0, CREST_X1,
            alpha=0.10, color="gray"
        )
        ax.set_xlim(DAM_X0, DAM_X1)
        ax.set_xlabel("Longitudinal distance x (m)")
        ax.set_ylabel("Bed elevation (m)")
        ax.set_title("B2 - Bed evolution profiles (0-600 s)")
        ax.legend()
        ax.grid(alpha=0.25)

        fig.tight_layout()
        fig.savefig(
            OUT / "01_bed_evolution_profiles_600s.png"
        )
        plt.close(fig)

        save_csv(
            OUT / "01_bed_evolution_profiles_600s.csv",
            ["time_s", "x_m", "y_m", "bed_elevation_m",
             "bed_change_m"],
            profile_rows
        )

        # Figure 2: cumulative erosion maps
        map_times = [60, 120, 300, 600]
        max_erosion = max(
            float(np.nanmax(
                np.maximum(-snapshots[t]["change"], 0.0)
            ))
            for t in map_times
        )

        levels = np.linspace(
            0.0, max(max_erosion, 1e-6), 21
        )

        fig, axes = plt.subplots(
            2, 2, figsize=(13, 8), sharex=True, sharey=True
        )

        for ax, target in zip(axes.flat, map_times):
            dz = snapshots[target]["change"]
            erosion = np.maximum(-dz, 0.0)

            cs = ax.tricontourf(
                triang, erosion,
                levels=levels,
                cmap="viridis",
                extend="max"
            )

            ax.axvline(
                CREST_X0, linestyle="--",
                color="cyan", linewidth=0.8
            )
            ax.axvline(
                CREST_X1, linestyle="--",
                color="cyan", linewidth=0.8
            )

            ax.set_xlim(DAM_X0 - 0.5, DAM_X1 + 0.5)
            ax.set_ylim(float(y.min()), float(y.max()))
            ax.set_title(
                f"t = {target} s | max erosion = "
                f"{np.nanmax(erosion):.4f} m"
            )
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
            ax.set_aspect("auto")

        fig.colorbar(
            cs,
            ax=axes.ravel().tolist(),
            label="Cumulative erosion depth (m)",
            shrink=0.85
        )

        fig.suptitle(
            "B2 - Cumulative vertical erosion evolution",
            fontsize=14
        )
        fig.subplots_adjust(
            left=0.07, right=0.89,
            bottom=0.08, top=0.89,
            hspace=0.30, wspace=0.18
        )
        fig.savefig(
            OUT / "02_bed_change_maps_600s.png"
        )
        plt.close(fig)

        # Figure 3: maximum erosion time series
        max_rows = []
        crest_rows = []

        for it, t in enumerate(times):
            dz = get_change(it)
            erosion = np.maximum(-dz, 0.0)

            max_rows.append([
                float(t),
                float(np.nanmax(erosion[dam_mask])),
                float(np.nanmax(dz[dam_mask]))
            ])

            crest_z = bed0[crest_mask] + dz[crest_mask]

            crest_rows.append([
                float(t),
                float(np.nanmean(crest_z)),
                float(np.nanmin(crest_z)),
                float(np.nanmax(crest_z))
            ])

        max_arr = np.asarray(max_rows, dtype=float)

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(
            max_arr[:, 0],
            max_arr[:, 1] * 1000,
            "-o", markersize=3
        )
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Maximum cumulative erosion (mm)")
        ax.set_title("B2 - Maximum erosion evolution")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(
            OUT / "03_maximum_erosion_600s.png"
        )
        plt.close(fig)

        save_csv(
            OUT / "03_maximum_erosion_600s.csv",
            ["time_s", "max_erosion_m", "max_deposition_m"],
            max_rows
        )

        # Figure 4: crest elevation
        crest_arr = np.asarray(crest_rows, dtype=float)

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(
            crest_arr[:, 0],
            crest_arr[:, 1],
            "-o", markersize=3,
            label="Mean crest elevation"
        )
        ax.plot(
            crest_arr[:, 0],
            crest_arr[:, 2],
            "--", linewidth=1,
            label="Minimum crest elevation"
        )

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Bed elevation (m)")
        ax.set_title("B2 - Crest elevation evolution")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            OUT / "04_crest_elevation_600s.png"
        )
        plt.close(fig)

        save_csv(
            OUT / "04_crest_elevation_600s.csv",
            ["time_s", "mean_crest_m",
             "min_crest_m", "max_crest_m"],
            crest_rows
        )

        # Summary
        print("")
        print("B2 600s ANALYSIS COMPLETE")
        print("")
        for target in TIMES:
            dz = snapshots[target]["change"]
            print(
                f"t={target:6.1f} s | "
                f"max erosion={np.max(-dz):.6f} m | "
                f"max deposition={np.max(dz):.6f} m"
            )

        print("")
        print("All PNG and CSV files saved in:")
        print(OUT)

    finally:
        slf.close()

if __name__ == "__main__":
    main()
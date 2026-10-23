# -*- coding: utf-8 -*-
"""
B2 - Tingsanchali & Chinnarasri (2001), Run B-2
Create the corrected full-flume 2D TELEMAC mesh and initial bed elevation.

Purpose
-------
This script builds the full experimental flume domain for the B2
erosion-process benchmark.

Experimental geometry used:
    Flume length            = 35.00 m
    Flume width             = 1.00 m
    Dam height              = 0.80 m
    Crest width             = 0.30 m
    Upstream slope          = 1V:3H
    Downstream slope        = 1V:5H

Local coordinate system:
    x = 0.00 m   : upstream/front edge of dam crest
    x = -17.50 m : upstream physical flume boundary
    x = -2.40 m  : upstream dam toe
    x = 0.30 m   : downstream edge of dam crest
    x = 4.30 m   : downstream dam toe
    x = 17.50 m  : downstream physical flume boundary

The full 35 m experimental flume is represented.
The upstream/front edge of the dam crest is taken as x = 0.00 m.

Mesh:
    dx = dy = 0.05 m
    structured points split into triangles

Expected mesh size:
    NX             = 701
    NY             = 21
    NPOIN          = 14721
    NELEM          = 28000
    Boundary nodes = 1440

IMPORTANT:
    This is B2 geometry construction only.
    No erosion is activated here.
    No TELEMAC simulation is run here.
"""

import os
import sys
import numpy as np


# ============================================================
# 1. PATHS
# ============================================================

ROOT = (
    r"F:\Telemac_Projects\Tuanzho_New"
    r"\16_B2_Tingsanchali_erosion_validation"
)

OUT_GEO_DIR = os.path.join(
    ROOT,
    "04_input",
)

OUT_MESH_DIR = os.path.join(
    ROOT,
    "02_mesh",
)

os.makedirs(
    OUT_GEO_DIR,
    exist_ok=True,
)

os.makedirs(
    OUT_MESH_DIR,
    exist_ok=True,
)

OUT_SLF = os.path.join(
    OUT_GEO_DIR,
    "B2_Tingsanchali_geo.slf",
)

OUT_NPZ = os.path.join(
    OUT_MESH_DIR,
    "B2_Tingsanchali_mesh_debug.npz",
)


# ============================================================
# 2. IMPORT TELEMAC PYTHON
# ============================================================

try:
    from data_manip.formats.selafin import Selafin

except Exception as exc:

    print("")
    print(
        "ERROR: TELEMAC Python modules could not be imported."
    )
    print(
        "Run this script from the configured "
        "openTELEMAC 9.1 PowerShell."
    )
    print("")
    print(exc)

    sys.exit(1)


# ============================================================
# 3. EXPERIMENTAL / COMPUTATIONAL PARAMETERS
# ============================================================

# ------------------------------------------------------------
# Full physical experimental flume
# ------------------------------------------------------------

XMIN = -17.50
XMAX = 17.50

YMIN = 0.00
YMAX = 1.00


# ------------------------------------------------------------
# Uniform mesh spacing
# ------------------------------------------------------------

DX = 0.05
DY = 0.05


# ------------------------------------------------------------
# Dam geometry
# ------------------------------------------------------------

DAM_HEIGHT = 0.80
CREST_WIDTH = 0.30

UPSTREAM_SLOPE_H = 3.0
DOWNSTREAM_SLOPE_H = 5.0


# ------------------------------------------------------------
# Derived x coordinates
# ------------------------------------------------------------

X_UPSTREAM_TOE = (
    -UPSTREAM_SLOPE_H
    * DAM_HEIGHT
)
# -2.40 m

X_CREST_FRONT = 0.00

X_CREST_BACK = CREST_WIDTH
# 0.30 m

X_DOWNSTREAM_TOE = (
    X_CREST_BACK
    + DOWNSTREAM_SLOPE_H
    * DAM_HEIGHT
)
# 4.30 m


# ============================================================
# 4. BASIC CONSISTENCY CHECKS
# ============================================================

def is_grid_aligned(
    value,
    origin,
    spacing,
    tol=1.0e-10,
):
    """
    Check whether a coordinate lies exactly on the
    structured computational grid.
    """

    q = (
        value
        - origin
    ) / spacing

    return (
        abs(
            q
            - round(q)
        )
        < tol
    )


important_x = [
    XMIN,
    X_UPSTREAM_TOE,
    X_CREST_FRONT,
    X_CREST_BACK,
    X_DOWNSTREAM_TOE,
    XMAX,
]


for xx in important_x:

    if not is_grid_aligned(
        xx,
        XMIN,
        DX,
    ):

        raise ValueError(
            f"Important x-coordinate {xx:.6f} m "
            f"is not aligned with DX={DX:.6f} m."
        )


if not is_grid_aligned(
    YMAX,
    YMIN,
    DY,
):

    raise ValueError(
        "YMAX is not aligned with DY."
    )


# Full-flume consistency checks

if abs(
    (XMAX - XMIN)
    - 35.00
) > 1.0e-12:

    raise ValueError(
        "The computational flume length is not 35.00 m."
    )


if abs(
    (YMAX - YMIN)
    - 1.00
) > 1.0e-12:

    raise ValueError(
        "The computational flume width is not 1.00 m."
    )


# ============================================================
# 5. CREATE GRID POINTS
# ============================================================

NX = (
    int(
        round(
            (XMAX - XMIN)
            / DX
        )
    )
    + 1
)

NY = (
    int(
        round(
            (YMAX - YMIN)
            / DY
        )
    )
    + 1
)


x_values = np.linspace(
    XMIN,
    XMAX,
    NX,
)

y_values = np.linspace(
    YMIN,
    YMAX,
    NY,
)


NPOIN = NX * NY


meshx = np.zeros(
    NPOIN,
    dtype=np.float64,
)

meshy = np.zeros(
    NPOIN,
    dtype=np.float64,
)


def node_id(
    ix,
    iy,
):
    """
    Zero-based node number.
    x changes fastest.
    """

    return (
        iy * NX
        + ix
    )


for iy, yy in enumerate(
    y_values
):

    for ix, xx in enumerate(
        x_values
    ):

        ip = node_id(
            ix,
            iy,
        )

        meshx[ip] = xx
        meshy[ip] = yy


# ============================================================
# 6. CREATE TRIANGULAR CONNECTIVITY
# ============================================================

triangles = []


for iy in range(
    NY - 1
):

    for ix in range(
        NX - 1
    ):

        n00 = node_id(
            ix,
            iy,
        )

        n10 = node_id(
            ix + 1,
            iy,
        )

        n01 = node_id(
            ix,
            iy + 1,
        )

        n11 = node_id(
            ix + 1,
            iy + 1,
        )

        # Alternate diagonal direction to avoid a systematic
        # directional bias across the whole flume.

        if (
            (ix + iy) % 2
            == 0
        ):

            triangles.append(
                [
                    n00,
                    n10,
                    n11,
                ]
            )

            triangles.append(
                [
                    n00,
                    n11,
                    n01,
                ]
            )

        else:

            triangles.append(
                [
                    n00,
                    n10,
                    n01,
                ]
            )

            triangles.append(
                [
                    n10,
                    n11,
                    n01,
                ]
            )


ikle = np.asarray(
    triangles,
    dtype=np.int32,
)


NELEM = ikle.shape[0]


# ============================================================
# 7. DEFINE INITIAL BED ELEVATION
# ============================================================

def initial_bottom(
    x
):
    """
    Initial B-2 dam geometry.

    z = 0 outside the embankment.

    Upstream slope:
        from (-2.40, 0.00)
        to   ( 0.00, 0.80)

        slope = 1V:3H

    Crest:
        x = 0.00 to 0.30
        z = 0.80

    Downstream slope:
        from (0.30, 0.80)
        to   (4.30, 0.00)

        slope = 1V:5H
    """

    if x < X_UPSTREAM_TOE:

        return 0.0

    elif x < X_CREST_FRONT:

        return (
            (x - X_UPSTREAM_TOE)
            / UPSTREAM_SLOPE_H
        )

    elif x <= X_CREST_BACK:

        return DAM_HEIGHT

    elif x <= X_DOWNSTREAM_TOE:

        return (
            DAM_HEIGHT
            - (
                (x - X_CREST_BACK)
                / DOWNSTREAM_SLOPE_H
            )
        )

    else:

        return 0.0


bottom = np.asarray(
    [
        initial_bottom(xx)
        for xx in meshx
    ],
    dtype=np.float64,
)


# Remove tiny floating-point noise

bottom[
    np.abs(bottom)
    < 1.0e-12
] = 0.0


# ============================================================
# 8. BUILD IPOBO
# ============================================================

# IPOBO = 0 for internal nodes.
# Boundary nodes receive consecutive positive integers.
#
# Counter-clockwise ordering:
#
#   bottom: left -> right
#   right : bottom -> top
#   top   : right -> left
#   left  : top -> bottom
#
# Corner nodes are included only once.


ipobo = np.zeros(
    NPOIN,
    dtype=np.int32,
)


boundary_nodes = []


# ------------------------------------------------------------
# Bottom boundary
# ------------------------------------------------------------

for ix in range(
    NX
):

    boundary_nodes.append(
        node_id(
            ix,
            0,
        )
    )


# ------------------------------------------------------------
# Right boundary, excluding bottom corner
# ------------------------------------------------------------

for iy in range(
    1,
    NY,
):

    boundary_nodes.append(
        node_id(
            NX - 1,
            iy,
        )
    )


# ------------------------------------------------------------
# Top boundary, excluding right corner
# ------------------------------------------------------------

for ix in range(
    NX - 2,
    -1,
    -1,
):

    boundary_nodes.append(
        node_id(
            ix,
            NY - 1,
        )
    )


# ------------------------------------------------------------
# Left boundary, excluding top and bottom corners
# ------------------------------------------------------------

for iy in range(
    NY - 2,
    0,
    -1,
):

    boundary_nodes.append(
        node_id(
            0,
            iy,
        )
    )


for ib, ip in enumerate(
    boundary_nodes,
    start=1,
):

    ipobo[ip] = ib


NBOUNDARY = len(
    boundary_nodes
)


# ============================================================
# 9. GEOMETRY CHECKS
# ============================================================

expected_npoin = (
    NX
    * NY
)

expected_nelem = (
    2
    * (NX - 1)
    * (NY - 1)
)

expected_boundary = (
    2 * NX
    + 2 * NY
    - 4
)


assert (
    NPOIN
    == expected_npoin
)

assert (
    NELEM
    == expected_nelem
)

assert (
    NBOUNDARY
    == expected_boundary
)


# ------------------------------------------------------------
# Explicit expected values for corrected full-flume mesh
# ------------------------------------------------------------

assert NX == 701
assert NY == 21

assert NPOIN == 14721
assert NELEM == 28000

assert NBOUNDARY == 1440


# ------------------------------------------------------------
# Bottom range
# ------------------------------------------------------------

assert (
    abs(
        np.min(bottom)
        - 0.0
    )
    < 1.0e-12
)

assert (
    abs(
        np.max(bottom)
        - DAM_HEIGHT
    )
    < 1.0e-12
)


# ------------------------------------------------------------
# Check representative centreline elevations
# ------------------------------------------------------------

def nearest_node(
    x_target,
    y_target=0.50,
):

    dist2 = (
        (meshx - x_target) ** 2
        + (meshy - y_target) ** 2
    )

    return int(
        np.argmin(
            dist2
        )
    )


check_locations = {
    "upstream_boundary": (
        -17.50,
        0.00,
    ),
    "upstream_toe": (
        -2.40,
        0.00,
    ),
    "crest_front": (
        0.00,
        0.80,
    ),
    "crest_back": (
        0.30,
        0.80,
    ),
    "downstream_toe": (
        4.30,
        0.00,
    ),
    "downstream_boundary": (
        17.50,
        0.00,
    ),
}


for name, (
    xx,
    zz_expected,
) in check_locations.items():

    ip = nearest_node(
        xx
    )

    zz = bottom[ip]

    if (
        abs(
            zz
            - zz_expected
        )
        > 1.0e-10
    ):

        raise RuntimeError(
            f"Geometry check failed at {name}: "
            f"x={xx}, "
            f"Z={zz}, "
            f"expected={zz_expected}"
        )


# ============================================================
# 10. TRIANGLE QUALITY / AREA CHECK
# ============================================================

x0 = meshx[
    ikle[:, 0]
]

y0 = meshy[
    ikle[:, 0]
]

x1 = meshx[
    ikle[:, 1]
]

y1 = meshy[
    ikle[:, 1]
]

x2 = meshx[
    ikle[:, 2]
]

y2 = meshy[
    ikle[:, 2]
]


triangle_areas = (
    0.5
    * np.abs(
        (x1 - x0)
        * (y2 - y0)
        - (x2 - x0)
        * (y1 - y0)
    )
)


MIN_TRIANGLE_AREA = float(
    np.min(
        triangle_areas
    )
)

MAX_TRIANGLE_AREA = float(
    np.max(
        triangle_areas
    )
)


N_ZERO_AREA = int(
    np.count_nonzero(
        triangle_areas
        <= 1.0e-14
    )
)


if N_ZERO_AREA != 0:

    raise RuntimeError(
        "Degenerate triangle(s) detected: "
        f"{N_ZERO_AREA}"
    )


# ============================================================
# 11. WRITE DEBUG NPZ
# ============================================================

np.savez(
    OUT_NPZ,
    meshx=meshx,
    meshy=meshy,
    ikle=ikle,
    ipobo=ipobo,
    bottom=bottom,
    boundary_nodes=np.asarray(
        boundary_nodes,
        dtype=np.int32,
    ),
    NX=NX,
    NY=NY,
    DX=DX,
    DY=DY,
    XMIN=XMIN,
    XMAX=XMAX,
    YMIN=YMIN,
    YMAX=YMAX,
    X_UPSTREAM_TOE=X_UPSTREAM_TOE,
    X_CREST_FRONT=X_CREST_FRONT,
    X_CREST_BACK=X_CREST_BACK,
    X_DOWNSTREAM_TOE=X_DOWNSTREAM_TOE,
)


# ============================================================
# 12. WRITE TELEMAC SELAFIN GEOMETRY
#
# IMPORTANT:
# This section retains the Selafin writing method used by the
# previous B2 geometry script that successfully produced the
# TELEMAC geometry file.
# ============================================================

slf = Selafin("")


slf.title = (
    "B2 TINGSANCHALI RUN B2 FULL FLUME GEOMETRY"
)


slf.nbv1 = 1
slf.nbv2 = 0


slf.varnames = [
    "BOTTOM          "
]

slf.varunits = [
    "M               "
]


slf.iparam = np.zeros(
    10,
    dtype=np.int32,
)


slf.nelem3 = NELEM
slf.npoin3 = NPOIN
slf.ndp3 = 3
slf.nplan = 1


slf.nelem2 = NELEM
slf.npoin2 = NPOIN
slf.ndp2 = 3


slf.ikle3 = ikle.copy()
slf.ikle2 = ikle.copy()


slf.ipob3 = ipobo.copy()
slf.ipob2 = ipobo.copy()


slf.meshx = meshx.copy()
slf.meshy = meshy.copy()


# Force double precision

slf.fformat = "SERAFIND"


slf.fole = {
    "hook": open(
        OUT_SLF,
        "wb",
    ),
    "name": OUT_SLF,
    "endian": ">",
    "float": (
        "d",
        8,
    ),
}


slf.append_header_slf()


# One geometry record at t = 0

slf.append_core_time_slf(
    0.0
)


values = np.zeros(
    (
        1,
        NPOIN,
    ),
    dtype=np.float64,
)

values[
    0,
    :
] = bottom


slf.append_core_vars_slf(
    values
)


slf.fole[
    "hook"
].close()


# ============================================================
# 13. REPORT
# ============================================================

print("")

print(
    "=" * 72
)

print(
    "B2 TINGSANCHALI CORRECTED FULL-FLUME "
    "GEOMETRY / MESH CREATED"
)

print(
    "=" * 72
)


print(
    f"Domain X              : "
    f"{XMIN:.2f} to {XMAX:.2f} m"
)

print(
    f"Domain Y              : "
    f"{YMIN:.2f} to {YMAX:.2f} m"
)

print(
    f"Total flume length    : "
    f"{XMAX - XMIN:.2f} m"
)

print(
    f"Flume width           : "
    f"{YMAX - YMIN:.2f} m"
)

print(
    f"DX, DY                : "
    f"{DX:.3f}, {DY:.3f} m"
)


print("")

print(
    f"NX                     : "
    f"{NX}"
)

print(
    f"NY                     : "
    f"{NY}"
)

print(
    f"NPOIN                  : "
    f"{NPOIN}"
)

print(
    f"NELEM                  : "
    f"{NELEM}"
)

print(
    f"Boundary nodes         : "
    f"{NBOUNDARY}"
)


print("")

print(
    "Dam geometry:"
)

print(
    f"  Upstream toe         : "
    f"x = {X_UPSTREAM_TOE:.2f} m"
)

print(
    f"  Crest front          : "
    f"x = {X_CREST_FRONT:.2f} m"
)

print(
    f"  Crest back           : "
    f"x = {X_CREST_BACK:.2f} m"
)

print(
    f"  Downstream toe       : "
    f"x = {X_DOWNSTREAM_TOE:.2f} m"
)

print(
    f"  Dam height           : "
    f"{DAM_HEIGHT:.2f} m"
)

print(
    f"  Crest width          : "
    f"{CREST_WIDTH:.2f} m"
)


print("")

print(
    f"BOTTOM min             : "
    f"{bottom.min():.6f} m"
)

print(
    f"BOTTOM max             : "
    f"{bottom.max():.6f} m"
)


print("")

print(
    "Centreline geometry checks:"
)


for name, (
    xx,
    zz_expected,
) in check_locations.items():

    ip = nearest_node(
        xx
    )

    print(
        f"  {name:20s}: "
        f"x={meshx[ip]:8.3f}, "
        f"y={meshy[ip]:6.3f}, "
        f"Z={bottom[ip]:8.5f}, "
        f"expected={zz_expected:8.5f}"
    )


print("")

print(
    "Triangle checks:"
)

print(
    f"  Minimum area         : "
    f"{MIN_TRIANGLE_AREA:.10f} m2"
)

print(
    f"  Maximum area         : "
    f"{MAX_TRIANGLE_AREA:.10f} m2"
)

print(
    f"  Zero-area triangles  : "
    f"{N_ZERO_AREA}"
)


print("")

print(
    "Files written:"
)

print(
    f"  SLF : {OUT_SLF}"
)

print(
    f"  NPZ : {OUT_NPZ}"
)


print("")

print(
    "STATUS: GEOMETRY / MESH CONSTRUCTION PASSED"
)

print(
    "=" * 72
)
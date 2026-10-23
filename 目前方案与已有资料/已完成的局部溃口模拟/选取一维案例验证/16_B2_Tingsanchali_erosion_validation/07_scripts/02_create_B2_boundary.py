# -*- coding: utf-8 -*-
"""
B2 - Tingsanchali & Chinnarasri (2001), Run B-2
Create TELEMAC boundary-condition file for the corrected full-flume mesh.

Geometry:
    x = -17.50 to 17.50 m
    y =  0.00 to  1.00 m

Boundary types:
    Left boundary (x = -17.50 m):
        prescribed discharge / liquid boundary
        LIHBOR = 4
        LIUBOR = 5
        LIVBOR = 5

    Right boundary (x = 17.50 m):
        free outflow / liquid boundary
        LIHBOR = 4
        LIUBOR = 4
        LIVBOR = 4

    Bottom and top boundaries:
        solid walls
        LIHBOR = 2
        LIUBOR = 2
        LIVBOR = 2

Corner treatment:
    Left two corners  -> inlet
    Right two corners -> outlet

IMPORTANT:
    This script only creates the CLI file.
    It does not run TELEMAC.
    It does not activate erosion.
"""

import os
import numpy as np


# ============================================================
# 1. PATHS
# ============================================================

ROOT = (
    r"F:\Telemac_Projects\Tuanzho_New"
    r"\16_B2_Tingsanchali_erosion_validation"
)

NPZ_FILE = os.path.join(
    ROOT,
    "02_mesh",
    "B2_Tingsanchali_mesh_debug.npz",
)

OUT_CLI = os.path.join(
    ROOT,
    "04_input",
    "B2_Tingsanchali.cli",
)


# ============================================================
# 2. DOMAIN PARAMETERS
# ============================================================

XMIN = -17.50
XMAX = 17.50

YMIN = 0.00
YMAX = 1.00

TOL = 1.0e-10


# ============================================================
# 3. TELEMAC BOUNDARY CODES
# ============================================================

# ------------------------------------------------------------
# Solid wall
# ------------------------------------------------------------

WALL_LIHBOR = 2
WALL_LIUBOR = 2
WALL_LIVBOR = 2


# ------------------------------------------------------------
# Upstream prescribed-discharge liquid boundary
#
# IMPORTANT:
# LIHBOR = 4 : free water depth / elevation
# LIUBOR = 5 : prescribed liquid-boundary velocity treatment
# LIVBOR = 5 : prescribed liquid-boundary velocity treatment
#
# This gives the required TELEMAC inlet code:
#
#                    4 5 5
#
# The actual Run B-2 inflow:
#
#       Q = 1.58 L/s = 0.00158 m3/s
#
# will be prescribed later in the TELEMAC steering file.
# ------------------------------------------------------------

INLET_LIHBOR = 4
INLET_LIUBOR = 5
INLET_LIVBOR = 5


# ------------------------------------------------------------
# Downstream free liquid boundary
# ------------------------------------------------------------

OUTLET_LIHBOR = 4
OUTLET_LIUBOR = 4
OUTLET_LIVBOR = 4


# ============================================================
# 4. LOAD VALIDATED MESH
# ============================================================

if not os.path.isfile(NPZ_FILE):

    raise FileNotFoundError(
        f"Mesh debug file not found:\n{NPZ_FILE}"
    )


data = np.load(NPZ_FILE)


meshx = data["meshx"]
meshy = data["meshy"]
ipobo = data["ipobo"]

boundary_nodes = data[
    "boundary_nodes"
].astype(np.int32)


NPOIN = meshx.size
NBOUNDARY = boundary_nodes.size


# ============================================================
# 5. BASIC CHECKS
# ============================================================

if NPOIN != 14721:

    raise RuntimeError(
        f"Unexpected NPOIN: {NPOIN}; "
        "expected 14721."
    )


if NBOUNDARY != 1440:

    raise RuntimeError(
        f"Unexpected boundary-node count: "
        f"{NBOUNDARY}; expected 1440."
    )


if np.count_nonzero(ipobo) != NBOUNDARY:

    raise RuntimeError(
        "IPOBO nonzero-node count does not match "
        "boundary-node count."
    )


# ------------------------------------------------------------
# Check that boundary_nodes follow IPOBO order 1...NBOUNDARY
# ------------------------------------------------------------

for ib, ip in enumerate(
    boundary_nodes,
    start=1,
):

    if int(ipobo[ip]) != ib:

        raise RuntimeError(
            "Boundary-node ordering does not match IPOBO: "
            f"boundary position={ib}, "
            f"node={ip + 1}, "
            f"IPOBO={int(ipobo[ip])}"
        )


# ============================================================
# 6. CLASSIFY BOUNDARY NODES
# ============================================================

records = []

count_inlet = 0
count_outlet = 0
count_wall = 0


for ib, ip in enumerate(
    boundary_nodes,
    start=1,
):

    x = float(meshx[ip])
    y = float(meshy[ip])

    # --------------------------------------------------------
    # Classification priority:
    #
    #   1. left liquid boundary
    #   2. right liquid boundary
    #   3. solid side walls
    #
    # Therefore:
    #   left corner nodes  -> inlet
    #   right corner nodes -> outlet
    # --------------------------------------------------------

    if abs(x - XMIN) <= TOL:

        boundary_type = "INLET"

        lihbor = INLET_LIHBOR
        liubor = INLET_LIUBOR
        livbor = INLET_LIVBOR

        count_inlet += 1


    elif abs(x - XMAX) <= TOL:

        boundary_type = "OUTLET"

        lihbor = OUTLET_LIHBOR
        liubor = OUTLET_LIUBOR
        livbor = OUTLET_LIVBOR

        count_outlet += 1


    elif (
        abs(y - YMIN) <= TOL
        or abs(y - YMAX) <= TOL
    ):

        boundary_type = "WALL"

        lihbor = WALL_LIHBOR
        liubor = WALL_LIUBOR
        livbor = WALL_LIVBOR

        count_wall += 1


    else:

        raise RuntimeError(
            "Boundary node could not be classified: "
            f"node={ip + 1}, "
            f"x={x:.12f}, "
            f"y={y:.12f}"
        )


    records.append(
        {
            "boundary_index": ib,
            "node_index": int(ip),
            "node_number": int(ip) + 1,
            "x": x,
            "y": y,
            "type": boundary_type,
            "LIHBOR": lihbor,
            "LIUBOR": liubor,
            "LIVBOR": livbor,
        }
    )


# ============================================================
# 7. EXPECTED BOUNDARY COUNTS
# ============================================================

# At each vertical liquid boundary:
#
# y = 0.00 ... 1.00 m
# DY = 0.05 m
#
# therefore:
#
# 1.00 / 0.05 + 1 = 21 nodes

EXPECTED_INLET = 21
EXPECTED_OUTLET = 21


# Remaining boundary nodes belong to the two side walls.
# Four corner nodes are already assigned to inlet/outlet.

EXPECTED_WALL = (
    NBOUNDARY
    - EXPECTED_INLET
    - EXPECTED_OUTLET
)


if count_inlet != EXPECTED_INLET:

    raise RuntimeError(
        f"Unexpected inlet-node count: "
        f"{count_inlet}; "
        f"expected {EXPECTED_INLET}."
    )


if count_outlet != EXPECTED_OUTLET:

    raise RuntimeError(
        f"Unexpected outlet-node count: "
        f"{count_outlet}; "
        f"expected {EXPECTED_OUTLET}."
    )


if count_wall != EXPECTED_WALL:

    raise RuntimeError(
        f"Unexpected wall-node count: "
        f"{count_wall}; "
        f"expected {EXPECTED_WALL}."
    )


if (
    count_inlet
    + count_outlet
    + count_wall
    != NBOUNDARY
):

    raise RuntimeError(
        "Boundary classification does not cover "
        "all boundary nodes."
    )


# ============================================================
# 8. WRITE TELEMAC CLI
# ============================================================

# TELEMAC boundary file columns:
#
# LIHBOR LIUBOR LIVBOR
# HBOR   UBOR   VBOR
# AUBOR  LITBOR TBOR
# ATBOR  BTBOR
# N      K
#
# Here:
#
# HBOR   = 0.0
# UBOR   = 0.0
# VBOR   = 0.0
# AUBOR  = 0.0
# LITBOR = 0
# TBOR   = 0.0
# ATBOR  = 0.0
# BTBOR  = 0.0
#
# N = global mesh node number (1-based)
# K = boundary-node number (1-based)


with open(
    OUT_CLI,
    "w",
    encoding="ascii",
    newline="\n",
) as f:

    for record in records:

        LIHBOR = record["LIHBOR"]
        LIUBOR = record["LIUBOR"]
        LIVBOR = record["LIVBOR"]

        HBOR = 0.0
        UBOR = 0.0
        VBOR = 0.0

        AUBOR = 0.0

        LITBOR = 0
        TBOR = 0.0

        ATBOR = 0.0
        BTBOR = 0.0

        N = record["node_number"]
        K = record["boundary_index"]

        line = (
            f"{LIHBOR:2d} "
            f"{LIUBOR:2d} "
            f"{LIVBOR:2d} "
            f"{HBOR:16.8f} "
            f"{UBOR:16.8f} "
            f"{VBOR:16.8f} "
            f"{AUBOR:16.8f} "
            f"{LITBOR:2d} "
            f"{TBOR:16.8f} "
            f"{ATBOR:16.8f} "
            f"{BTBOR:16.8f} "
            f"{N:8d} "
            f"{K:8d}\n"
        )

        f.write(line)


# ============================================================
# 9. VERIFY WRITTEN CLI
# ============================================================

with open(
    OUT_CLI,
    "r",
    encoding="ascii",
) as f:

    written_lines = [
        line
        for line in f
        if line.strip()
    ]


if len(written_lines) != NBOUNDARY:

    raise RuntimeError(
        f"Written CLI line count = "
        f"{len(written_lines)}, "
        f"expected {NBOUNDARY}."
    )


verify_inlet = 0
verify_outlet = 0
verify_wall = 0


for line in written_lines:

    parts = line.split()


    if len(parts) != 13:

        raise RuntimeError(
            "Unexpected number of CLI columns: "
            f"{len(parts)}"
        )


    lihbor = int(parts[0])
    liubor = int(parts[1])
    livbor = int(parts[2])

    node_number = int(parts[11])
    boundary_index = int(parts[12])


    if (
        boundary_index < 1
        or boundary_index > NBOUNDARY
    ):

        raise RuntimeError(
            "Invalid boundary index in written CLI."
        )


    ip = node_number - 1

    x = float(meshx[ip])
    y = float(meshy[ip])


    # --------------------------------------------------------
    # Verify inlet = 4 5 5
    # --------------------------------------------------------

    if abs(x - XMIN) <= TOL:

        if (
            lihbor,
            liubor,
            livbor,
        ) != (
            4,
            5,
            5,
        ):

            raise RuntimeError(
                "Incorrect inlet boundary code "
                "in written CLI. "
                f"Node={node_number}, "
                f"found=({lihbor},{liubor},{livbor}), "
                "expected=(4,5,5)."
            )

        verify_inlet += 1


    # --------------------------------------------------------
    # Verify outlet = 4 4 4
    # --------------------------------------------------------

    elif abs(x - XMAX) <= TOL:

        if (
            lihbor,
            liubor,
            livbor,
        ) != (
            4,
            4,
            4,
        ):

            raise RuntimeError(
                "Incorrect outlet boundary code "
                "in written CLI. "
                f"Node={node_number}, "
                f"found=({lihbor},{liubor},{livbor}), "
                "expected=(4,4,4)."
            )

        verify_outlet += 1


    # --------------------------------------------------------
    # Verify walls = 2 2 2
    # --------------------------------------------------------

    elif (
        abs(y - YMIN) <= TOL
        or abs(y - YMAX) <= TOL
    ):

        if (
            lihbor,
            liubor,
            livbor,
        ) != (
            2,
            2,
            2,
        ):

            raise RuntimeError(
                "Incorrect wall boundary code "
                "in written CLI. "
                f"Node={node_number}, "
                f"found=({lihbor},{liubor},{livbor}), "
                "expected=(2,2,2)."
            )

        verify_wall += 1


    else:

        raise RuntimeError(
            "Written CLI contains a node "
            "that is not on the domain boundary."
        )


if verify_inlet != count_inlet:

    raise RuntimeError(
        "Written CLI inlet count mismatch."
    )


if verify_outlet != count_outlet:

    raise RuntimeError(
        "Written CLI outlet count mismatch."
    )


if verify_wall != count_wall:

    raise RuntimeError(
        "Written CLI wall count mismatch."
    )


# ============================================================
# 10. FINAL CODE CHECKS
# ============================================================

if (
    INLET_LIHBOR,
    INLET_LIUBOR,
    INLET_LIVBOR,
) != (
    4,
    5,
    5,
):

    raise RuntimeError(
        "INLET code is not 4 5 5."
    )


if (
    OUTLET_LIHBOR,
    OUTLET_LIUBOR,
    OUTLET_LIVBOR,
) != (
    4,
    4,
    4,
):

    raise RuntimeError(
        "OUTLET code is not 4 4 4."
    )


if (
    WALL_LIHBOR,
    WALL_LIUBOR,
    WALL_LIVBOR,
) != (
    2,
    2,
    2,
):

    raise RuntimeError(
        "WALL code is not 2 2 2."
    )


# ============================================================
# 11. REPORT
# ============================================================

print("")

print("=" * 72)

print(
    "B2 TINGSANCHALI CORRECTED FULL-FLUME "
    "BOUNDARY CREATED"
)

print("=" * 72)


print(
    f"Mesh nodes             : "
    f"{NPOIN}"
)

print(
    f"Boundary nodes         : "
    f"{NBOUNDARY}"
)


print("")

print(
    "Boundary classification:"
)

print(
    f"  Inlet x={XMIN:6.2f} m   : "
    f"{count_inlet} nodes "
    f"(4 5 5)"
)

print(
    f"  Outlet x={XMAX:6.2f} m  : "
    f"{count_outlet} nodes "
    f"(4 4 4)"
)

print(
    f"  Solid side walls     : "
    f"{count_wall} nodes "
    f"(2 2 2)"
)


print("")

print(
    "Corner-node treatment:"
)

print(
    "  Left two corners     : INLET (4 5 5)"
)

print(
    "  Right two corners    : OUTLET (4 4 4)"
)


print("")

print(
    "CLI verification:"
)

print(
    f"  Written lines        : "
    f"{len(written_lines)}"
)

print(
    f"  Verified inlet       : "
    f"{verify_inlet}"
)

print(
    f"  Verified outlet      : "
    f"{verify_outlet}"
)

print(
    f"  Verified walls       : "
    f"{verify_wall}"
)


print("")

print(
    "Boundary codes:"
)

print(
    "  INLET                : 4 5 5"
)

print(
    "  OUTLET               : 4 4 4"
)

print(
    "  WALL                 : 2 2 2"
)


print("")

print(
    "Boundary file:"
)

print(
    f"  {OUT_CLI}"
)


print("")

print(
    "STATUS: BOUNDARY / CLI CONSTRUCTION PASSED"
)

print("=" * 72)
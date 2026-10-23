import os
import re
import sys

# ============================================================
# B2 Tingsanchali tau_c sensitivity-case generator
#
# Purpose:
#   Generate two additional tau_c sensitivity cases from the
#   already validated erosion reference run.
#
# Fixed parameters:
#   p     = 0.40
#   k_d   = 5.0e-4 m/(Pa s)
#
# tau_c cases:
#   0.001 Pa  -> generated
#   0.010 Pa  -> existing reference run; NOT regenerated
#   0.100 Pa  -> generated
#
# IMPORTANT:
#   This script only creates Fortran and CAS input files.
#   It does NOT run TELEMAC.
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

REF_FORTRAN = os.path.join(
    RUN_DIR,
    "user_condin_h_erosion_ref.f"
)

REF_CAS = os.path.join(
    RUN_DIR,
    "B2_Tingsanchali_erosion_ref.cas"
)


# ------------------------------------------------------------
# CASE DEFINITIONS
# ------------------------------------------------------------

CASES = [
    {
        "tag": "tauc_0p001",
        "tau_c": 0.001,
        "fortran_value": "0.001D0",
    },
    {
        "tag": "tauc_0p100",
        "tau_c": 0.100,
        "fortran_value": "0.100D0",
    },
]


# ------------------------------------------------------------
# CHECK REFERENCE FILES
# ------------------------------------------------------------

print()
print("=" * 72)
print("B2 TAU_C SENSITIVITY CASE GENERATOR")
print("=" * 72)

if not os.path.isfile(REF_FORTRAN):
    print("ERROR: Reference Fortran file not found:")
    print(REF_FORTRAN)
    sys.exit(1)

if not os.path.isfile(REF_CAS):
    print("ERROR: Reference CAS file not found:")
    print(REF_CAS)
    sys.exit(1)

print("Reference Fortran:")
print(" ", REF_FORTRAN)

print("Reference CAS:")
print(" ", REF_CAS)

print()


# ------------------------------------------------------------
# READ VALIDATED REFERENCE FILES
# ------------------------------------------------------------

with open(
    REF_FORTRAN,
    "r",
    encoding="utf-8"
) as f:
    ref_fortran_text = f.read()

with open(
    REF_CAS,
    "r",
    encoding="utf-8"
) as f:
    ref_cas_text = f.read()


# ------------------------------------------------------------
# SAFETY CHECK:
# Make sure the validated reference parameters are present.
# ------------------------------------------------------------

tau_pattern = re.compile(
    r"^\s*TAU_C\s*=\s*0\.01D0\s*$",
    flags=re.MULTILINE
)

kd_pattern = re.compile(
    r"^\s*K_D\s*=\s*5\.D-4\s*$",
    flags=re.MULTILINE
)

porosity_pattern = re.compile(
    r"^\s*POROSITY\s*=\s*0\.40D0\s*$",
    flags=re.MULTILINE
)


tau_matches = tau_pattern.findall(
    ref_fortran_text
)

kd_matches = kd_pattern.findall(
    ref_fortran_text
)

porosity_matches = porosity_pattern.findall(
    ref_fortran_text
)


if len(tau_matches) != 1:
    print(
        "ERROR: Expected exactly one active "
        "TAU_C = 0.01D0 line."
    )
    print(
        "Found:",
        len(tau_matches)
    )
    sys.exit(1)

if len(kd_matches) != 1:
    print(
        "ERROR: Expected exactly one active "
        "K_D = 5.D-4 line."
    )
    print(
        "Found:",
        len(kd_matches)
    )
    sys.exit(1)

if len(porosity_matches) != 1:
    print(
        "ERROR: Expected exactly one active "
        "POROSITY = 0.40D0 line."
    )
    print(
        "Found:",
        len(porosity_matches)
    )
    sys.exit(1)


print("Reference parameter check:")
print("  POROSITY = 0.40")
print("  TAU_C    = 0.01 Pa")
print("  K_D      = 5.0e-4 m/(Pa s)")
print("  STATUS   = OK")
print()


# ------------------------------------------------------------
# GENERATE CASES
# ------------------------------------------------------------

generated_files = []


for case in CASES:

    tag = case["tag"]
    tau_c = case["tau_c"]
    tau_value = case["fortran_value"]


    # --------------------------------------------------------
    # OUTPUT NAMES
    # --------------------------------------------------------

    fortran_name = (
        f"user_condin_h_erosion_{tag}.f"
    )

    cas_name = (
        f"B2_Tingsanchali_erosion_{tag}.cas"
    )

    result_name = (
        f"B2_Tingsanchali_erosion_{tag}_res.slf"
    )


    fortran_path = os.path.join(
        RUN_DIR,
        fortran_name
    )

    cas_path = os.path.join(
        RUN_DIR,
        cas_name
    )


    # --------------------------------------------------------
    # SAFETY:
    # Do not silently overwrite an existing generated case.
    # --------------------------------------------------------

    if os.path.exists(fortran_path):
        print("ERROR: File already exists:")
        print(fortran_path)
        print("No file was overwritten.")
        sys.exit(1)

    if os.path.exists(cas_path):
        print("ERROR: File already exists:")
        print(cas_path)
        print("No file was overwritten.")
        sys.exit(1)


    # --------------------------------------------------------
    # MODIFY FORTRAN
    #
    # ONLY ACTIVE TAU_C ASSIGNMENT IS CHANGED.
    # K_D, porosity, limiter, wet-count logic, erosion zone,
    # initialization and TELEMAC interface remain identical.
    # --------------------------------------------------------

    new_fortran_text, n_tau = (
        tau_pattern.subn(
            f"      TAU_C = {tau_value}",
            ref_fortran_text
        )
    )

    if n_tau != 1:
        print(
            "ERROR: TAU_C replacement count != 1"
        )
        sys.exit(1)


    # Update descriptive text so the generated file is
    # self-identifying. This does not affect calculations.

    new_fortran_text = new_fortran_text.replace(
        "REFERENCE RUN - NOT FINAL CALIBRATION",
        (
            "TAU_C SENSITIVITY RUN - "
            "NOT FINAL CALIBRATION"
        )
    )


    # --------------------------------------------------------
    # WRITE FORTRAN
    # --------------------------------------------------------

    with open(
        fortran_path,
        "w",
        encoding="utf-8",
        newline="\n"
    ) as f:

        f.write(
            new_fortran_text
        )


    # --------------------------------------------------------
    # MODIFY CAS
    #
    # Change only:
    #   TITLE
    #   RESULTS FILE
    #   FORTRAN FILE
    #
    # All validated hydrodynamic settings remain unchanged.
    # --------------------------------------------------------

    new_cas_text = ref_cas_text


    new_cas_text, n_title = re.subn(
        r"^TITLE\s*=.*$",
        (
            "TITLE = "
            f"'B2 TINGSANCHALI TAUC {tau_c:.3f} PA'"
        ),
        new_cas_text,
        count=1,
        flags=re.MULTILINE
    )


    new_cas_text, n_result = re.subn(
        r"^RESULTS FILE\s*=.*$",
        (
            "RESULTS FILE = "
            f"'{result_name}'"
        ),
        new_cas_text,
        count=1,
        flags=re.MULTILINE
    )


    new_cas_text, n_fortran = re.subn(
        r"^FORTRAN FILE\s*=.*$",
        (
            "FORTRAN FILE = "
            f"'{fortran_name}'"
        ),
        new_cas_text,
        count=1,
        flags=re.MULTILINE
    )


    if n_title != 1:
        print(
            "ERROR: TITLE replacement count != 1"
        )
        sys.exit(1)

    if n_result != 1:
        print(
            "ERROR: RESULTS FILE replacement count != 1"
        )
        sys.exit(1)

    if n_fortran != 1:
        print(
            "ERROR: FORTRAN FILE replacement count != 1"
        )
        sys.exit(1)


    # --------------------------------------------------------
    # WRITE CAS
    # --------------------------------------------------------

    with open(
        cas_path,
        "w",
        encoding="utf-8",
        newline="\n"
    ) as f:

        f.write(
            new_cas_text
        )


    generated_files.append(
        (
            tau_c,
            fortran_name,
            cas_name,
            result_name
        )
    )


# ------------------------------------------------------------
# FINAL REPORT
# ------------------------------------------------------------

print("=" * 72)
print("GENERATED TAU_C CASES")
print("=" * 72)

for (
    tau_c,
    fortran_name,
    cas_name,
    result_name
) in generated_files:

    print()
    print(
        f"tau_c = {tau_c:.3f} Pa"
    )

    print(
        "  Fortran :",
        fortran_name
    )

    print(
        "  CAS     :",
        cas_name
    )

    print(
        "  Result  :",
        result_name
    )


print()
print("-" * 72)

print(
    "Existing reference case:"
)

print(
    "  tau_c = 0.010 Pa"
)

print(
    "  B2_Tingsanchali_erosion_ref.cas"
)

print(
    "  B2_Tingsanchali_erosion_ref_res.slf"
)

print("-" * 72)

print()
print(
    "IMPORTANT: TELEMAC has NOT been run by this script."
)

print(
    "STATUS: TAU_C CASE GENERATION PASSED"
)

print("=" * 72)
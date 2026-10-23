import os
import re
import sys


# ============================================================
# B2 Tingsanchali - k_d calibration/sensitivity case generator
#
# Validated reference case:
#   porosity = 0.40
#   tau_c    = 0.01 Pa
#   k_d      = 5.0e-4 m/(Pa s)
#
# Target k_d set:
#   0.5e-4  -> generated
#   1.0e-4  -> generated
#   2.5e-4  -> generated
#   5.0e-4  -> existing reference; NOT regenerated
#   8.0e-4  -> generated
#
# IMPORTANT:
#   1. The validated reference files are READ ONLY.
#   2. Only K_D and identifying file names are changed.
#   3. tau_c, porosity, hydrodynamics, limiter, wet-count
#      logic, erosion zone and initial conditions remain fixed.
#   4. This script does NOT run TELEMAC.
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
#
# Naming:
#
#   kd_0p5e4 -> 0.5e-4 m/(Pa s)
#   kd_1p0e4 -> 1.0e-4 m/(Pa s)
#   kd_2p5e4 -> 2.5e-4 m/(Pa s)
#   kd_8p0e4 -> 8.0e-4 m/(Pa s)
#
# Reference:
#
#   kd = 5.0e-4 m/(Pa s)
#   already exists and is NOT generated again.
# ------------------------------------------------------------

CASES = [
    {
        "tag": "kd_0p5e4",
        "kd": 0.5e-4,
        "fortran_value": "0.5D-4",
    },
    {
        "tag": "kd_1p0e4",
        "kd": 1.0e-4,
        "fortran_value": "1.0D-4",
    },
    {
        "tag": "kd_2p5e4",
        "kd": 2.5e-4,
        "fortran_value": "2.5D-4",
    },
    {
        "tag": "kd_8p0e4",
        "kd": 8.0e-4,
        "fortran_value": "8.0D-4",
    },
]


# ------------------------------------------------------------
# HEADER
# ------------------------------------------------------------

print()
print("=" * 78)
print("B2 K_D CASE GENERATOR")
print("=" * 78)


# ------------------------------------------------------------
# CHECK REFERENCE FILES
# ------------------------------------------------------------

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
# STRICT REFERENCE-PARAMETER CHECK
# ------------------------------------------------------------
#
# We require the mother file to contain exactly:
#
#   POROSITY = 0.40D0
#   TAU_C    = 0.01D0
#   K_D      = 5.D-4
#
# If the mother file has changed unexpectedly, STOP.
# ------------------------------------------------------------

porosity_pattern = re.compile(
    r"^\s*POROSITY\s*=\s*0\.40D0\s*$",
    flags=re.MULTILINE
)

tau_pattern = re.compile(
    r"^\s*TAU_C\s*=\s*0\.01D0\s*$",
    flags=re.MULTILINE
)

kd_pattern = re.compile(
    r"^\s*K_D\s*=\s*5\.D-4\s*$",
    flags=re.MULTILINE
)


porosity_matches = porosity_pattern.findall(
    ref_fortran_text
)

tau_matches = tau_pattern.findall(
    ref_fortran_text
)

kd_matches = kd_pattern.findall(
    ref_fortran_text
)


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


print("Reference parameter check:")
print("  POROSITY = 0.40")
print("  TAU_C    = 0.01 Pa")
print("  K_D      = 5.0e-4 m/(Pa s)")
print("  STATUS   = OK")
print()


# ------------------------------------------------------------
# CHECK CAS REFERENCE NAMES
# ------------------------------------------------------------

required_cas_strings = [
    "user_condin_h_erosion_ref.f",
    "B2_Tingsanchali_erosion_ref_res.slf",
]


for required_text in required_cas_strings:

    count = ref_cas_text.count(
        required_text
    )

    if count != 1:
        print(
            "ERROR: Expected exactly one occurrence in "
            "reference CAS:"
        )
        print(
            " ",
            required_text
        )
        print(
            "Found:",
            count
        )
        sys.exit(1)


print("Reference CAS filename check:")
print("  STATUS   = OK")
print()


# ------------------------------------------------------------
# GENERATE CASES
# ------------------------------------------------------------

generated_files = []


for case in CASES:

    tag = case["tag"]
    kd = case["kd"]
    kd_value = case["fortran_value"]


    # --------------------------------------------------------
    # UNIQUE OUTPUT NAMES
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
    # DO NOT OVERWRITE EXISTING FILES
    # --------------------------------------------------------

    if os.path.exists(fortran_path):
        print()
        print("ERROR: Generated Fortran file already exists:")
        print(fortran_path)
        print()
        print("No existing file was overwritten.")
        sys.exit(1)


    if os.path.exists(cas_path):
        print()
        print("ERROR: Generated CAS file already exists:")
        print(cas_path)
        print()
        print("No existing file was overwritten.")
        sys.exit(1)


    # --------------------------------------------------------
    # MODIFY FORTRAN
    #
    # ONLY:
    #   K_D = ...
    #
    # is numerically changed.
    # --------------------------------------------------------

    new_fortran_text, n_kd = (
        kd_pattern.subn(
            f"      K_D = {kd_value}",
            ref_fortran_text
        )
    )


    if n_kd != 1:
        print(
            "ERROR: K_D replacement count != 1"
        )
        sys.exit(1)


    # Update descriptive comment only.
    # No numerical effect.

    new_fortran_text = new_fortran_text.replace(
        "REFERENCE RUN - NOT FINAL CALIBRATION",
        (
            "K_D CALIBRATION RUN - "
            "NOT FINAL CALIBRATION"
        )
    )


    # --------------------------------------------------------
    # POST-MODIFICATION PARAMETER SAFETY CHECK
    # --------------------------------------------------------

    if len(
        porosity_pattern.findall(
            new_fortran_text
        )
    ) != 1:
        print(
            "ERROR: POROSITY changed unexpectedly."
        )
        sys.exit(1)


    if len(
        tau_pattern.findall(
            new_fortran_text
        )
    ) != 1:
        print(
            "ERROR: TAU_C changed unexpectedly."
        )
        sys.exit(1)


    expected_kd_pattern = re.compile(
        r"^\s*K_D\s*=\s*"
        + re.escape(kd_value)
        + r"\s*$",
        flags=re.MULTILINE
    )


    if len(
        expected_kd_pattern.findall(
            new_fortran_text
        )
    ) != 1:
        print(
            "ERROR: Generated K_D value check failed:"
        )
        print(
            " ",
            kd_value
        )
        sys.exit(1)


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
    # Only:
    #   TITLE
    #   RESULTS FILE
    #   FORTRAN FILE
    #
    # All numerical TELEMAC settings stay unchanged.
    # --------------------------------------------------------

    new_cas_text = ref_cas_text


    new_cas_text, n_title = re.subn(
        r"^TITLE\s*=.*$",
        (
            "TITLE = "
            f"'B2 TINGSANCHALI KD {kd:.1e}'"
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
    # VERIFY CAS
    # --------------------------------------------------------

    if new_cas_text.count(
        fortran_name
    ) != 1:
        print(
            "ERROR: Generated FORTRAN FILE check failed."
        )
        sys.exit(1)


    if new_cas_text.count(
        result_name
    ) != 1:
        print(
            "ERROR: Generated RESULTS FILE check failed."
        )
        sys.exit(1)


    if (
        "user_condin_h_erosion_ref.f"
        in new_cas_text
    ):
        print(
            "ERROR: Generated CAS still references "
            "reference Fortran."
        )
        sys.exit(1)


    if (
        "B2_Tingsanchali_erosion_ref_res.slf"
        in new_cas_text
    ):
        print(
            "ERROR: Generated CAS still references "
            "reference result file."
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
        {
            "kd": kd,
            "fortran": fortran_name,
            "cas": cas_name,
            "result": result_name,
        }
    )


# ------------------------------------------------------------
# FINAL REPORT
# ------------------------------------------------------------

print("=" * 78)
print("GENERATED K_D CASES")
print("=" * 78)


for item in generated_files:

    print()

    print(
        f"k_d = {item['kd']:.1e} m/(Pa s)"
    )

    print(
        "  Fortran :",
        item["fortran"]
    )

    print(
        "  CAS     :",
        item["cas"]
    )

    print(
        "  Result  :",
        item["result"]
    )


print()
print("-" * 78)

print(
    "Existing reference case:"
)

print(
    "  k_d   = 5.0e-4 m/(Pa s)"
)

print(
    "  tau_c = 0.01 Pa"
)

print(
    "  B2_Tingsanchali_erosion_ref.cas"
)

print(
    "  B2_Tingsanchali_erosion_ref_res.slf"
)

print("-" * 78)

print()

print(
    "IMPORTANT:"
)

print(
    "  Four new input cases were generated."
)

print(
    "  The validated reference case was not modified."
)

print(
    "  TELEMAC has NOT been run."
)

print()

print(
    "STATUS: B2 K_D CASE GENERATION PASSED"
)

print("=" * 78)
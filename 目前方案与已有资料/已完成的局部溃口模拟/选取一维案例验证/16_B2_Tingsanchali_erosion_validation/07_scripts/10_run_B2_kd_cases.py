import os
import sys
import subprocess
from datetime import datetime


# ============================================================
# B2 Tingsanchali - automated k_d TELEMAC runner
#
# Runs ONLY the four newly generated k_d cases:
#
#   0.5e-4
#   1.0e-4
#   2.5e-4
#   8.0e-4  m/(Pa s)
#
# Existing reference:
#
#   5.0e-4 m/(Pa s)
#
# is NOT rerun.
#
# Behaviour:
#   1. Check all required CAS files.
#   2. Run cases sequentially through PowerShell.
#   3. Save complete stdout/stderr for every case.
#   4. Require:
#        return code = 0
#        "CORRECT END OF RUN" in log
#        expected SLF result exists
#   5. STOP immediately if any case fails.
#
# IMPORTANT:
#   Run this script from an openTELEMAC-configured
#   PowerShell environment.
#
# Windows/openTELEMAC note:
#   telemac2d.py is normally resolved correctly by the
#   configured PowerShell environment, but direct
#   subprocess execution with shell=False may fail to
#   resolve the .py launcher.
#
#   Therefore this runner explicitly invokes:
#
#       powershell.exe -NoProfile -Command
#           telemac2d.py <case.cas>
#
#   The current configured environment is inherited by
#   that PowerShell process.
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

OUTPUT_DIR = os.path.join(
    ROOT,
    "06_output"
)

LOG_DIR = os.path.join(
    OUTPUT_DIR,
    "kd_run_logs"
)

os.makedirs(
    LOG_DIR,
    exist_ok=True
)


# ------------------------------------------------------------
# CASES
# ------------------------------------------------------------

CASES = [
    {
        "label": "KD050",
        "kd": 0.5e-4,
        "cas": "B2_Tingsanchali_erosion_kd_0p5e4.cas",
        "result": "B2_Tingsanchali_erosion_kd_0p5e4_res.slf",
    },
    {
        "label": "KD100",
        "kd": 1.0e-4,
        "cas": "B2_Tingsanchali_erosion_kd_1p0e4.cas",
        "result": "B2_Tingsanchali_erosion_kd_1p0e4_res.slf",
    },
    {
        "label": "KD250",
        "kd": 2.5e-4,
        "cas": "B2_Tingsanchali_erosion_kd_2p5e4.cas",
        "result": "B2_Tingsanchali_erosion_kd_2p5e4_res.slf",
    },
    {
        "label": "KD800",
        "kd": 8.0e-4,
        "cas": "B2_Tingsanchali_erosion_kd_8p0e4.cas",
        "result": "B2_Tingsanchali_erosion_kd_8p0e4_res.slf",
    },
]


# ------------------------------------------------------------
# HEADER
# ------------------------------------------------------------

print()
print("=" * 78)
print("B2 AUTOMATED K_D TELEMAC RUNNER")
print("=" * 78)

print("Run directory:")
print(" ", RUN_DIR)

print("Log directory:")
print(" ", LOG_DIR)

print()


# ------------------------------------------------------------
# PRE-FLIGHT CHECK
# ------------------------------------------------------------

print("PRE-FLIGHT CHECK")
print("-" * 78)


for case in CASES:

    cas_path = os.path.join(
        RUN_DIR,
        case["cas"]
    )

    if not os.path.isfile(cas_path):

        print("ERROR: CAS file not found:")
        print(" ", cas_path)

        sys.exit(1)

    print(
        f"{case['label']:6s} | "
        f"k_d = {case['kd']:.1e} | "
        f"CAS OK"
    )


print("-" * 78)
print("All four CAS files found.")
print()


# ------------------------------------------------------------
# CHECK POWERSHELL + TELEMAC COMMAND
# ------------------------------------------------------------

print("TELEMAC COMMAND CHECK")
print("-" * 78)


try:

    command_check = subprocess.run(
        [
            "powershell.exe",
            "-NoProfile",
            "-Command",
            "Get-Command telemac2d.py | "
            "Select-Object -ExpandProperty Source"
        ],
        cwd=RUN_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace"
    )

except FileNotFoundError:

    print(
        "ERROR: powershell.exe was not found."
    )

    sys.exit(1)


telemac_command_path = (
    command_check.stdout.strip()
)


if (
    command_check.returncode != 0
    or not telemac_command_path
):

    print(
        "ERROR: telemac2d.py cannot be resolved "
        "inside the inherited PowerShell environment."
    )

    print()
    print("PowerShell output:")
    print(command_check.stdout)

    sys.exit(1)


print(
    "telemac2d.py resolved as:"
)

print(
    " ",
    telemac_command_path
)

print(
    "STATUS: TELEMAC command available"
)

print("-" * 78)
print()


# ------------------------------------------------------------
# RUN CASES
# ------------------------------------------------------------

completed = []


for number, case in enumerate(
    CASES,
    start=1
):

    label = case["label"]
    kd = case["kd"]

    cas_name = case["cas"]
    result_name = case["result"]

    result_path = os.path.join(
        RUN_DIR,
        result_name
    )

    log_name = (
        f"{label}_run.log"
    )

    log_path = os.path.join(
        LOG_DIR,
        log_name
    )


    print()
    print("=" * 78)

    print(
        f"CASE {number}/{len(CASES)} : "
        f"{label}"
    )

    print(
        f"k_d = {kd:.1e} m/(Pa s)"
    )

    print(
        f"CAS = {cas_name}"
    )

    print("=" * 78)


    # --------------------------------------------------------
    # SAFETY:
    # Do not silently overwrite an existing result.
    # --------------------------------------------------------

    if os.path.isfile(result_path):

        print()
        print(
            "ERROR: Expected result file already exists:"
        )

        print(
            " ",
            result_path
        )

        print()
        print(
            "The runner will not overwrite an existing "
            "simulation result."
        )

        print(
            "Resolve this file manually before rerunning."
        )

        sys.exit(1)


    # --------------------------------------------------------
    # RUN TELEMAC THROUGH POWERSHELL
    # --------------------------------------------------------

    start_time = datetime.now()

    print(
        "Start:",
        start_time.strftime(
            "%Y-%m-%d %H:%M:%S"
        )
    )

    print(
        "Running TELEMAC..."
    )


    ps_command = (
        f"telemac2d.py '.\\{cas_name}'"
    )


    try:

        process = subprocess.run(
            [
                "powershell.exe",
                "-NoProfile",
                "-Command",
                ps_command
            ],
            cwd=RUN_DIR,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace"
        )

    except Exception as exc:

        print()
        print(
            "ERROR while launching TELEMAC:"
        )

        print(
            str(exc)
        )

        sys.exit(1)


    end_time = datetime.now()

    elapsed = (
        end_time
        - start_time
    ).total_seconds()


    # --------------------------------------------------------
    # SAVE COMPLETE LOG
    # --------------------------------------------------------

    with open(
        log_path,
        "w",
        encoding="utf-8",
        errors="replace"
    ) as f:

        f.write(
            process.stdout
        )


    print(
        "End:",
        end_time.strftime(
            "%Y-%m-%d %H:%M:%S"
        )
    )

    print(
        f"Elapsed = {elapsed:.1f} s"
    )

    print(
        "Return code =",
        process.returncode
    )

    print(
        "Log saved:"
    )

    print(
        " ",
        log_path
    )


    # --------------------------------------------------------
    # RUN SUCCESS CHECKS
    # --------------------------------------------------------

    correct_end = (
        "CORRECT END OF RUN"
        in process.stdout
    )

    result_exists = os.path.isfile(
        result_path
    )


    print()
    print("POST-RUN CHECK")

    print(
        "  Return code       :",
        process.returncode
    )

    print(
        "  CORRECT END OF RUN:",
        correct_end
    )

    print(
        "  Result SLF exists :",
        result_exists
    )


    # --------------------------------------------------------
    # FAILURE -> STOP THE WHOLE BATCH
    # --------------------------------------------------------

    if (
        process.returncode != 0
        or not correct_end
        or not result_exists
    ):

        print()
        print("=" * 78)

        print(
            f"FAILED: {label}"
        )

        print(
            "The batch has been stopped."
        )

        print(
            "No later k_d cases will be run."
        )

        print()
        print(
            "Inspect log:"
        )

        print(
            " ",
            log_path
        )

        print()
        print(
            "Last TELEMAC output:"
        )

        print("-" * 78)

        output_lines = (
            process.stdout.splitlines()
        )

        for line in output_lines[-40:]:
            print(line)

        print("-" * 78)

        sys.exit(1)


    # --------------------------------------------------------
    # EXTRACT FINAL EROSION DIAGNOSTIC FROM LOG
    # --------------------------------------------------------

    diagnostic_lines = []

    capture = False

    for line in process.stdout.splitlines():

        if (
            "B-2 EROSION CHECK AT ABOUT 30 s"
            in line
        ):
            capture = True

        if capture:

            diagnostic_lines.append(
                line
            )

        if (
            capture
            and "SUM DZ/STEP"
            in line
        ):
            break


    print()
    print("FINAL 30-S EROSION DIAGNOSTIC")
    print("-" * 78)

    if diagnostic_lines:

        for line in diagnostic_lines:
            print(line)

    else:

        print(
            "WARNING: 30-s erosion diagnostic "
            "block not found in log."
        )


    # --------------------------------------------------------
    # RECORD SUCCESS
    # --------------------------------------------------------

    completed.append(
        label
    )

    print()
    print(
        f"STATUS: {label} PASSED"
    )


# ------------------------------------------------------------
# FINAL BATCH REPORT
# ------------------------------------------------------------

print()
print("=" * 78)
print("B2 K_D BATCH RUN COMPLETED")
print("=" * 78)

print(
    "Completed new cases:"
)

for label in completed:

    print(
        " ",
        label
    )


print()
print(
    "Existing reference case was NOT rerun:"
)

print(
    "  KD500 = 5.0e-4 m/(Pa s)"
)

print()

print(
    "Complete k_d set:"
)

print(
    "  KD050 = 0.5e-4"
)

print(
    "  KD100 = 1.0e-4"
)

print(
    "  KD250 = 2.5e-4"
)

print(
    "  KD500 = 5.0e-4  [existing reference]"
)

print(
    "  KD800 = 8.0e-4"
)

print()

print(
    "STATUS: ALL NEW B2 K_D CASES PASSED"
)

print("=" * 78)
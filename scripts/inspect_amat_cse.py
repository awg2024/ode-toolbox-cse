import json
import os
from pprint import pprint

import odetoolbox


# ============================================================
# 1. FIND THE REPOSITORY ROOT
# ============================================================

# This script lives inside:
#
# ode-toolbox-cse/scripts/inspect_amat_cse.py
#
# So:
# dirname(__file__)                  -> scripts/
# dirname(dirname(__file__))         -> ode-toolbox-cse/
#
SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

REPO_ROOT = os.path.dirname(
    SCRIPT_DIR
)


# ============================================================
# 2. DEFINE INPUT / OUTPUT LOCATIONS
# ============================================================

INPUT_FILE = os.path.join(
    REPO_ROOT,
    "scripts",
    "tests",
    "amat_ode_toolbox_input.json",
)

OUTPUT_DIR = os.path.join(
    REPO_ROOT,
    "outputs",
)


# Create outputs/ if it does not exist.
os.makedirs(
    OUTPUT_DIR,
    exist_ok=True,
)


BASELINE_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "amat_solver_baseline.json",
)

CSE_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "amat_solver_cse.json",
)

CSE_CONDITIONS_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "amat_solver_cse_conditions.json",
)


# ============================================================
# 3. LOAD THE AMAT ODE-TOOLBOX INPUT
# ============================================================

print("\nLoading AMAT ODE-toolbox input from:")
print(INPUT_FILE)


with open(
    INPUT_FILE,
    "r",
    encoding="utf-8",
) as file:

    model = json.load(file)


print("\nAMAT input loaded successfully.")


# ============================================================
# 4. RUN BASELINE
# ============================================================

print(
    "\n"
    "============================================================\n"
    "RUN 1: BASELINE — NO CSE\n"
    "============================================================"
)


baseline = odetoolbox.analysis(
    model,

    # We do not need PyGSL for this experiment.
    disable_stiffness_check=True,

    # Ordinary ODE-toolbox behaviour.
    enable_cse=False,
)


# ============================================================
# 5. RUN CSE WITHOUT CONDITION-BRANCH OPTIMISATION
# ============================================================

print(
    "\n"
    "============================================================\n"
    "RUN 2: CSE — CONDITIONS LEFT UNTOUCHED\n"
    "============================================================"
)


cse = odetoolbox.analysis(
    model,

    disable_stiffness_check=True,

    enable_cse=True,

    # Important:
    # singularity branches remain exactly as generated
    # by ODE-toolbox.
    enable_cse_condition_branches=False,
)


# ============================================================
# 6. RUN CSE INCLUDING INDEPENDENT CONDITION BRANCHES
# ============================================================

print(
    "\n"
    "============================================================\n"
    "RUN 3: CSE — CONDITION BRANCH CSE ENABLED\n"
    "============================================================"
)


cse_conditions = odetoolbox.analysis(
    model,

    disable_stiffness_check=True,

    enable_cse=True,

    # Each singularity branch gets its OWN local CSE.
    enable_cse_condition_branches=True,
)


# ============================================================
# 7. HELPER FUNCTION FOR SAVING JSON
# ============================================================

def save_json(data, filename):

    with open(
        filename,
        "w",
        encoding="utf-8",
    ) as file:

        json.dump(
            data,
            file,
            indent=2,
            sort_keys=True,
        )

    print(
        "Saved:",
        filename,
    )


# ============================================================
# 8. SAVE ALL THREE OUTPUTS
# ============================================================

print(
    "\n"
    "============================================================\n"
    "SAVING OUTPUT FILES\n"
    "============================================================"
)


save_json(
    baseline,
    BASELINE_OUTPUT,
)


save_json(
    cse,
    CSE_OUTPUT,
)


save_json(
    cse_conditions,
    CSE_CONDITIONS_OUTPUT,
)


# ============================================================
# 9. HELPER FOR PRINTING SOLVER STRUCTURE
# ============================================================

def inspect_solver_collection(
    title,
    solver_collection,
):

    print(
        "\n"
        "============================================================"
    )

    print(title)

    print(
        "============================================================"
    )

    for index, solver in enumerate(
        solver_collection
    ):

        print(
            f"\nSolver {index}"
        )

        print(
            "Solver type:",
            solver.get("solver"),
        )

        print(
            "State variables:",
            solver.get(
                "state_variables",
                [],
            ),
        )

        print(
            "Top-level keys:",
            list(
                solver.keys()
            ),
        )

        # ----------------------------------------------------
        # Top-level CSE
        # ----------------------------------------------------

        if "cse" in solver:

            print(
                "\nTop-level CSE replacements:"
            )

            pprint(
                solver["cse"]
            )

        else:

            print(
                "\nNo top-level CSE replacements."
            )

        # ----------------------------------------------------
        # Conditional analytical solver
        # ----------------------------------------------------

        if "conditions" in solver:

            print(
                "\nSingularity conditions:"
            )

            for condition, branch in (
                solver[
                    "conditions"
                ].items()
            ):

                print(
                    f"\n  Condition: {condition}"
                )

                print(
                    "  Branch keys:",
                    list(
                        branch.keys()
                    ),
                )

                if "cse" in branch:

                    print(
                        "  Branch-local CSE:"
                    )

                    pprint(
                        branch["cse"]
                    )

                else:

                    print(
                        "  No branch-local CSE."
                    )


# ============================================================
# 10. PRINT THE THREE RESULTS
# ============================================================

inspect_solver_collection(
    "BASELINE",
    baseline,
)


inspect_solver_collection(
    "CSE — CONDITIONS UNTOUCHED",
    cse,
)


inspect_solver_collection(
    "CSE — CONDITION BRANCHES ENABLED",
    cse_conditions,
)


# ============================================================
# 11. FINAL SUMMARY
# ============================================================

print(
    "\n"
    "============================================================\n"
    "DONE\n"
    "============================================================"
)

print(
    "\nBaseline JSON:"
)

print(
    BASELINE_OUTPUT
)


print(
    "\nCSE JSON:"
)

print(
    CSE_OUTPUT
)


print(
    "\nCSE + condition branches JSON:"
)

print(
    CSE_CONDITIONS_OUTPUT
)
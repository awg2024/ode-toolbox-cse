import json
import os
from pprint import pprint
import odetoolbox

# we need to adjust these paths for the diagnostics utils 
SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

REPO_ROOT = os.path.dirname(
    SCRIPT_DIR
)

# out and input dirs 
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

INPUT_FILE = r"json_testfile/amat_ode_toolbox_input.json"
CSE_OUPUT = r"outputs/"

print("\nLoading AMAT ODE-toolbox input from:")
print(INPUT_FILE)


with open(INPUT_FILE,"r",encoding="utf-8") as file: # read only, utf-8 reads text characters 

    model = json.load(file) # load json 

print("\nAMAT input loaded successfully.")


# baseline no cse implemented 
print("baseline amat building, no cse activated. ")

baseline = odetoolbox.analysis(
    model,
    # We do not need PyGSL for this experiment.
    disable_stiffness_check=True,
    # Ordinary ODE-toolbox behaviour.
    disable_singularity_mitigation=True,
    enable_cse=False,
)

# print("run 2 - cse, conditions left untouched")
# cse = odetoolbox.analysis(
#     model,
#     disable_stiffness_check=True,
#     enable_cse=True,
#     enable_cse_condition_branches=False, # don't touch conitions 
# )

# print("run 3 cse, singularity flag turned on")
# cse_conditions = odetoolbox.analysis(
#     model,
#     disable_stiffness_check=True,
#     enable_cse=True,
#     # Each singularity branch gets its OWN local CSE.
#     enable_cse_condition_branches=True,
# )

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
            sort_keys=True)

    print("Saved:",filename)


print("save output files")


save_json(baseline,BASELINE_OUTPUT)
save_json(cse,CSE_OUTPUT)
save_json(cse_conditions,CSE_CONDITIONS_OUTPUT)

def inspect_solver_collection(
    title,
    solver_collection):

    print(title)
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
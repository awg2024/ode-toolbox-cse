import json
import odetoolbox

path = r"/p/project1/paj2623/gray2/ode-toolbox-cse/scripts/tests/amat_ode_toolbox_input.json"

with open(
    path,
    "r",
) as file:

    model = json.load(file)


result = odetoolbox.analysis(
    model,
    disable_stiffness_check=True,
    enable_cse=False,
)

print("amat solver classification")

for index, solver in enumerate(result):

    print(f"\nSolver {index}")
    print("type:",solver["solver"])
    print("states:")

    for state in solver["state_variables"]:

        print("   ",state)

    print("conditions:" "conditions" in solver)
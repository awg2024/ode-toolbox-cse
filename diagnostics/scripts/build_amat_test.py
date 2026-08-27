import json
import os
from pprint import pprint
import odetoolbox

INPUT_FILE = r"json_testfile/amat_ode_toolbox_input.json"
CSE_OUPUT = r"outputs/"
CSE_CONDITIONS_OUTPUT = r"outputs/"

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
    disable_singularity_mitigation=True,  # disable searches for alternative solvers for singularities in the odetoolbox 
    enable_cse=False,
)

# disable_singularity_mitigation=True, # disable searches for singularities in the odetoolbox 

print("run 2 - cse, conditions left untouched")
cse = odetoolbox.analysis(
    model,
    disable_stiffness_check=True,
    enable_cse=True,
    disable_singularity_mitigation=True,  # disable searches for alternative solvers for singularities in the odetoolbox 
    enable_cse_condition_branches=False, # don't touch conitions 
)

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


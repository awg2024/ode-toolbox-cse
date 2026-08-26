"""
Run script: time python scripts/test_amat_combination_explosion.py

This evaluate the combinational explosion for the amat compilation 


"""

import json
import time 
import odetoolbox

path = r"/p/project1/paj2623/gray2/ode-toolbox-cse/scripts/tests/amat_ode_toolbox_input.json"

with open(
    path,
    "r",
) as file:

    model = json.load(file)

    print("generatirng one analytical amat solver instead of the 64 combinations")
    start = time.perf_counter()

    result = odetoolbox.analysis(model, disable_stiffness_check=True, disable_singularity_detection=True, enable_cse=False)

    elapsed = (time.perf_counter() - start)

    print(f"completed in: {elapsed}")

    for index, solver in enumerate(result):

        print(f"Solver: {index}")
        print("type:", solver["solver"])
        print("state:", solver["state_variables"])
        print("condiitons:", "conditions" in solver)


"""
this script is for inspecting condition statements from neuronal models and augementing the cse optimiser so we can correctly 
optimise the eq's without breaking them. 
"""

import sympy 
from odetoolbox.expression_optimisation import (common_subexpression_elimination, restore_cse_expression)
import odetoolbox

def print_solver_structure(
    solver,
    indent=4):

    prefix = " " * indent
    print(prefix, "solver:", solver.get("solver")) # print formatting indent. 

    print(prefix, "keys:", list(solver.keys())) # print keys and headers 

    if "propagators" in solver: # searches for propagator in high level dict 

        print(prefix, "propagators:",
            list(solver["propagators"].keys())) 

    if "update_expressions" in solver: # searches for update_expression in high level dict 

        print(prefix,"update expressions:",
            list(solver["update_expressions"].keys()))

    if "conditions" in solver: # searches for singularity condition 

        for condition, branch in (
            solver["conditions"].items()):

            print(prefix, f"condition: {condition}") # print key and items of condition

            print(prefix, " branch keys:", list(branch.keys()))

result = odetoolbox.analysis(
    'amat2_psc_exp',  # wont' we need an analytically solved nestml to test? 
    disable_stiffness_check=True,
    disable_singularity_detection=False,
    disable_singularity_mitigation=False,
    enable_cse=False,
)

for solver in result:
    print_solver_structure(
        solver
    )
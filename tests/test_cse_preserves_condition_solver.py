"""
In this test we will see how the cse solver reacts to a condition in the solver
"""

import sympy
from odetoolbox.expression_optimisation import (apply_cse_to_solver)
import json
import odetoolbox

def test_cse_preserves_condition_solver():

    indict = { # defining indict non-linear that will need a singularity condition 
        "dynamics": [
            {
                "expression": "V_m' = -V_m / tau_m + I/C_m",
                "initial_value": "0"
            }
        ],
        "parameters": {
            "tau_m": "10",
            "C_m": "250"
        },
        "options": {
            "output_timestep_symbol": "__h"
        }
    }


    solver = odetoolbox.analysis(indict)  # baseline, no CSE
    result = odetoolbox.analysis(indict, enable_cse=True) # apply cse 
    
    print(json.dumps(result, indent=2)) # print result output
    
    assert result[0]["conditions"] == solver[0]["conditions"] # ensure singularity conditions are present in both 

    propagators = result[0]["conditions"]["default"]["propagators"]
    update_exprs = result[0]["conditions"]["default"]["update_expressions"]
    all_text = " ".join(list(propagators.values()) + list(update_exprs.values()))
    assert "cse" not in all_text  # or whatever the actual generated symbol prefix is
    
    
 
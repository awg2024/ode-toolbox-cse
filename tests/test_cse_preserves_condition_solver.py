import sympy
from odetoolbox.expression_optimisation import (apply_cse_to_solver)

def test_cse_preserves_condition_solver():

    x, y = sympy.symbols("x y",real=True)

    solver = { # define solver as indict that would be passed into indict 
        "solver": "analytical",

        "state_variables": [
            "x",
            "y",
        ],

        "conditions": { #  apply_cse_to_solver only looks for high level dict levels (propagator, update) 

            "default": {
                "propagators": {
                    "P_x": sympy.exp(-x),
                    "P_y": sympy.exp(-y),
                },

                "update_expressions": {
                    "x": x,
                    "y": y,
                },
            },

            "(tau_1 == tau_2)": {
                "propagators": {
                    "P_x": 1 + x,
                    "P_y": 1 + y,
                },

                "update_expressions": {
                    "x": x,
                    "y": y,
                },
            },
        },
    }

    result = apply_cse_to_solver(solver) # apply cse to solver 

    assert result["conditions"] == solver["conditions"] # checking the result and solver are kept the same 

    assert "cse" not in result # make sure cse is not applied. 
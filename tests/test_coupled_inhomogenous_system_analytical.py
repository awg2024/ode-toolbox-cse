"""
This test is created after the PR #107 to test that we can correctly identify and compile analytical indicts in ODEToolBox
"""


from .context import odetoolbox

def test_coupled_inhomogeneous_system_is_analytical():

    model = {
        "dynamics": [
            {
            "expression":
                "x' = -x / tau_x + I", # if we set to initial values x(0) y(0), I still stands alone driving the eq (inhomogenous)
            "initial_value":
                "0",
            },
            {
            "expression":
                "y' = x - y / tau_y", # the value of the second eq requires the value of x to calculate y' (coupled)
            "initial_value":
                "0"
            },
        ],

        "parameters": {
            "tau_x": "10",
            "tau_y": "20",
            "I": "1",
        }
    }

    result = odetoolbox.analysis(model, disable_stiffness_check=True, enable_cse=False) # so we are expecting the odetoolbox to correctly classify this model as linear 
    
    assert len(result) == 1  
    solver = result[0]
    assert solver["solver"] == "analytical" # if analytical is part of the output we break 
    assert set(solver["state_variables"]) == {"x","y"} 
    
import json
import odetoolbox

def test_analysis_cse_is_json_serializable():

    indict = {
        "dynamics": [
            {
                "expression":
                    "x' = x + exp(x + y)",
                "initial_value":
                    "0",
            },
            {
                "expression":
                    "y' = y + 2 * exp(x + y)",
                "initial_value":
                    "0",
            },
        ]
    }

    result = odetoolbox.analysis(
        indict,
        disable_analytic_solver=True,
        disable_stiffness_check=True,
        enable_cse=True,
    )

    # This should raise nothing.
    encoded = json.dumps(result)

    assert isinstance(
        encoded,
        str,
    )
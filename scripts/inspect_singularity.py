import json
import odetoolbox

indict = {
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

result = odetoolbox.analysis(indict)   # or: disable_singularity_detection=True to compare
print(json.dumps(result, indent=2))
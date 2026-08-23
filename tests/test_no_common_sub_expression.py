
import sympy
from odetoolbox.expression_optimisation import (apply_cse_to_solver)

def test_cse_no_common_expression():

    x, y = sympy.symbols(
        "x y",
        real=True,
    )

    solver = {
        "solver": "numeric",
        "update_expressions": {
            "x": x + 1, #  optimised expressions 
            "y": y * 2,
        },
    }

    result = apply_cse_to_solver(solver)

    assert "cse" not in result #   check cse hasn't been conducted anyway 
    assert (result["update_expressions"] == solver["update_expressions"]) #  make sure update expressions are kept the same from result to solver 
import json
import sympy
from odetoolbox import (
    _read_global_config,
    _from_json_to_shapes,
    _find_analytically_solvable_equations,
)

from odetoolbox.system_of_shapes import (
    SystemOfShapes,
)


"""
This script is used as a utility helper definition area
to help build the amat system and test it as it passes through 
the odetoolbox. 
"""

def build_amat_system(
    input_file="json_testfile/amat_ode_toolbox_input.json",
):
    """
    Load the captured AMAT ODE-toolbox input and
    reproduce the system construction performed
    inside odetoolbox._analysis().

    No propagator is generated here.
    """

    with open(
        input_file,
        "r",
        encoding="utf-8",
    ) as file:
        model = json.load(file)

    # Apply ODE-toolbox configuration from JSON.
    _read_global_config(model)

    # Convert parameter names to real SymPy symbols,
    # matching _analysis().
    parameters = {}

    for name, value in model.get(
        "parameters",
        {},
    ).items():

        if isinstance(name, str):
            symbol = sympy.Symbol(
                name,
                real=True,
            )
        else:
            symbol = name

        parameters[symbol] = value

    # Parse dynamics into Shape objects.
    shapes, parameters = (
        _from_json_to_shapes(
            model,
            parameters=parameters,
        )
    )

    # Construct x' = A*x + b + c.
    shape_sys = (
        SystemOfShapes.from_shapes(
            shapes,
            parameters=parameters,
        )
    )

    # Apply the current analytical classification
    # algorithm, including PR #107.
    _, classification = (
        _find_analytically_solvable_equations(
            shape_sys,
            shapes,
            parameters=parameters,
        )
    )

    analytic_symbols = [
        symbol
        for symbol, is_analytic
        in classification.items()
        if is_analytic
    ]

    analytic_sys = (
        shape_sys.get_sub_system(
            analytic_symbols
        )
    )

    return (
        model,
        shape_sys,
        analytic_sys,
        classification,
    )
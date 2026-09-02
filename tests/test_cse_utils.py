#
# cse_test_utils.py
#
# This file is part of the NEST ODE toolbox.
#
# Copyright (C) 2017 The NEST Initiative
#
# The NEST ODE toolbox is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 2 of
# the License, or (at your option) any later version.
#
# The NEST ODE toolbox is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Helpers for CSE tests in isolated ODE-toolbox. 

"""

from __future__ import annotations

import copy
import json

from dataclasses import dataclass
from pathlib import Path

import sympy

from .context import odetoolbox
import os

@dataclass
class CSEAnalysisPair: 
    """
    Baseline and CSE-enabled results prdocued from the same ODE input
    """
    
    # defining output from conditions from _analysis
    baseline_solver: list
    cse_sovler: list

    baseline_shape_sys: list
    cse_shape_sys: list

    baseline_shapes: list
    cse_shape: list 


def load_test_json(filename: str) -> dict:
    """
    Load a JSON test fixture file by name.
    Resolved relative to this file's directory, not the CWD.
    """
    filepath = Path(__file__).parent / filename

    with open(filepath, "r") as f:
        return json.load(f)

def run_cse_analysis_pair(indict, **analysis_kwargs):

    """
    Run exactly the same ODE-toolbox twice (cse on/off)
    """

    # ensuring enable_cse has not been manipulated in keyword arguments 
    if "enable_cse" in analysis_kwargs: 
        raise ValueError("run_cse_analysis_pair controls enable_cse itself")

    baseline_solvers, baseline_shape_sys, baseline_shapes = odetoolbox._analysis(indict, enable_cse=False, **analysis_kwargs)

    cse_solvers, cse_shape_sys, cse_shapes = odetoolbox._analysis(indict, enable_cse=True, **analysis_kwargs)

    return CSEAnalysisPair(
        baseline_solvers=baseline_solvers,
        baseline_shape_sys=baseline_shape_sys,
        baseline_shapes=baseline_shapes,
        cse_solvers=cse_solvers,
        cse_shape_sys=cse_shape_sys,
        cse_shapes=cse_shapes)

def get_solver(solvers, solver_type):

    """
    retrieve exactly one solver by type (numeric, numeric-rk4, numeric-4kf45)
    """

    if solver_type == "numeric":
        matches = [solver for solver in solvers if solver["solver"].startswith("numeric")]
    else:
        matches = [solver for solver in solvers if sovler["solver"] == solver_type]

    assert len(matches) == 1, f"expected exactly one {solver_type} solver"

    return matches[0]


# Sympify can execute arbitary string inputs 
# pass in local dict to ensure expressions are correctly identified 
_SYMPY_LOCALS = {
    "exp": sympy.exp,
    "sin": sympy.sin,
    "cos": sympy.cos,
    "tan": sympy.tan,
    "log": sympy.log,
    "sqrt": sympy.sqrt,
    "Abs": sympy.Abs,
    "Min": sympy.Min,
    "Max": sympy.Max,
    "Piecewise": sympy.Piecewise}


def parse_expression(expression, extra_locals=None):
    """
    Convert serialized ODE-toolbox output back to a SymPy expression.
    """

    if isinstance(expression, sympy.Basic):
        return expression

    locals_dict = dict(_SYMPY_LOCALS)

    if extra_locals:
        locals_dict.update(extra_locals)

    return sympy.sympify(
        expression,
        locals=locals_dict)


def deserialize_replacements(replacements):
    
    """
    Convert serialized CSE metadata back into ordered SymPy replacement pairs.

    Supports the current format:

        {
            "__ode_cse_update_0": "x + y",
            "__ode_cse_update_1": "exp(__ode_cse_update_0)"
        }

    and also accepts raw ``[(Symbol, Expr), ...]`` data for convenience.
    """

    if replacements is None:
        return []

    #
    # Already raw SymPy replacements.
    #
    if isinstance(replacements, (list, tuple)):

        if all(
            isinstance(item, tuple)
            and len(item) == 2
            for item in replacements
        ):
            return list(replacements)

        #
        # Compatibility with possible
        # [{"symbol": ..., "expression": ...}] representations? this is the representation for codeg? 
        #
        if all(
            isinstance(item, dict)
            and "symbol" in item
            and "expression" in item
            for item in replacements
        ):
            replacements = {
                item["symbol"]: item["expression"]
                for item in replacements
            }

    if not isinstance(replacements, dict):
        raise TypeError(
            "Unknown CSE replacement representation: "
            f"{type(replacements)}"
        )

    #
    # Create every temporary Symbol first.
    #
    temporary_symbols = {
        symbol_name: sympy.Symbol(symbol_name)
        for symbol_name in replacements.keys()
    }

    result = []

    #
    # Dictionary insertion order preserves SymPy CSE dependency order.
    #
    for symbol_name, expression in replacements.items():

        temporary_symbol = temporary_symbols[symbol_name]

        parsed_expression = parse_expression(expression, extra_locals=temporary_symbols)

        result.append((temporary_symbol, parsed_expression))

    return result


def restore_cse_expression(
    reduced_expression,
    replacements,
):
    """
    Reconstruct a CSE-reduced expression by substituting every generated
    temporary variable back into the expression.

    Replacements are processed in reverse because later temporaries may depend
    on earlier temporaries. Ensures the mathematical expression has not been broken during CSE.
    """

    replacements = deserialize_replacements(replacements)

    temporary_locals = {
        str(symbol): symbol
        for symbol, _ in replacements
    }

    restored = parse_expression(
        reduced_expression,
        extra_locals=temporary_locals,
    )

    for temporary, expression in reversed(
        replacements
    ):
        restored = restored.subs(
            temporary,
            expression,
        )

    return sympy.simplify(restored)


def assert_expressions_equivalent(expected, actual):
    
    """
    Assert symbolic mathematical equivalence.
    """

    expected = parse_expression(expected)
    actual = parse_expression(actual)

    difference = sympy.simplify(expected - actual)

    if difference == 0:
        return # test has passed 

    equals_result = expected.equals(actual)

    assert equals_result is True, (
        "\nExpressions are not equivalent.\n"
        f"Expected:\n{expected}\n\n"
        f"Actual:\n{actual}\n\n"
        f"Difference:\n{difference}\n"
    )


def assert_cse_region_equivalent(baseline_solver, cse_solver, region_name, *, require_cse=True):
    """
    Prove mathematical equivalence between one baseline execution region and
    its CSE-reduced counterpart.

    Examples:
        region_name="propagators"
        region_name="update_expressions"
    """

    assert region_name in baseline_solver
    assert region_name in cse_solver

    baseline_region = baseline_solver[region_name]

    cse_region = cse_solver[region_name]

    assert set(baseline_region.keys()) == set(cse_region.keys())

    cse_metadata = (
        cse_solver
        .get("cse", {})
        .get(region_name)
    )

    if require_cse:
        assert cse_metadata, (
            f"Expected CSE replacements in "
            f"{region_name}, but none were generated"
        )

    for expression_name in baseline_region:

        baseline_expression = (
            baseline_region[expression_name]
        )

        cse_expression = (
            cse_region[expression_name]
        )

        if cse_metadata:

            reconstructed_expression = (
                restore_cse_expression(
                    cse_expression,
                    cse_metadata,
                )
            )

        else:
            reconstructed_expression = (
                parse_expression(cse_expression))
            
        assert_expressions_equivalent(baseline_expression, reconstructed_expression)


def assert_cse_region_profitable(
    baseline_solver,
    cse_solver,
    region_name,
):
    """
    Assert that the retained CSE transformation actually reduces symbolic
    operation count. 
    """

    replacements_serialized = (
        cse_solver
        .get("cse", {})
        .get(region_name)
    )

    assert replacements_serialized

    replacements = deserialize_replacements(replacements_serialized)

    temporary_locals = {
        str(symbol): symbol
        for symbol, _ in replacements
    }

    baseline_cost = sum(
        int(
            sympy.count_ops(
                parse_expression(expression)
            )
        )
        for expression in baseline_solver[
            region_name
        ].values()
    )

    replacement_cost = sum(
        int(sympy.count_ops(expression))
        for _, expression in replacements
    )

    reduced_cost = sum(
        int(
            sympy.count_ops(
                parse_expression(
                    expression,
                    extra_locals=temporary_locals,
                )
            )
        )
        for expression in cse_solver[
            region_name
        ].values()
    )

    if baseline_cost > 0:
        percent_reduction = ((baseline_cost - cse_cost) / baseline_cost) * 100
        percent_str = f"{percent_reduction:.2f}% reduction"
    else:
        percent_str = "0.00% change (baseline is 0)"

    print(f"[CSE INFO] Cost BEFORE: {baseline_cost}, AFTER: {cse_cost}, Change: {percent_str}")


    assert cse_cost < baseline_cost, (
        f"CSE was retained for {region_name} "
        f"but did not lower operation count: "
        f"{baseline_cost} -> {cse_cost}"
    )


def assert_solver_metadata_preserved(
    baseline_solver,
    cse_solver,
):
    """
    CSE is an expression optimisation only. It must not alter solver identity,
    state variables, parameters, or initial values.
    """

    metadata_keys = (
        "solver",
        "state_variables",
        "initial_values",
        "parameters",
    )

    for key in metadata_keys:

        if (
            key in baseline_solver
            or key in cse_solver
        ):
            assert (
                baseline_solver.get(key)
                == cse_solver.get(key)
            ), (
                f"CSE unexpectedly changed "
                f"solver metadata field: {key}"
            )


def assert_shape_structure_preserved(pair):
    """
    CSE happens after Shape/SystemOfShapes construction and therefore must not
    change the internal ODE representation.
    """

    baseline_symbols = [
        str(symbol)
        for symbol in pair.baseline_shape_sys.x_
    ]

    cse_symbols = [
        str(symbol)
        for symbol in pair.cse_shape_sys.x_
    ]

    assert baseline_symbols == cse_symbols

    baseline_shapes = [
        (
            str(shape.symbol),
            shape.order,
        )
        for shape in pair.baseline_shapes
    ]

    cse_shapes = [
        (
            str(shape.symbol),
            shape.order,
        )
        for shape in pair.cse_shapes
    ]

    assert baseline_shapes == cse_shapes


def assert_cse_serialized(
    solver,
):
    """
    Ensure the public/internal analysis result contains JSON-safe CSE metadata.
    """

    if "cse" in solver:

        for region_name, replacements in (
            solver["cse"].items()
        ):

            assert isinstance(
                replacements,
                dict,
            ), (
                f"CSE metadata for {region_name} "
                "was not serialized"
            )

            assert all(
                isinstance(symbol, str)
                for symbol in replacements.keys()
            )

            assert all(
                isinstance(expression, str)
                for expression
                in replacements.values()
            )

    for condition_solver in (
        solver.get(
            "conditions",
            {}
        ).values()
    ):
        assert_cse_serialized(
            condition_solver
        )


def _collect_cse_symbols(solver):
    """
    Collect all temporary CSE names in one solver block.
    """

    result = set()

    for replacements in (
        solver.get("cse", {}).values()
    ):

        if isinstance(replacements, dict):
            result.update(
                replacements.keys()
            )
        else:
            result.update(
                str(symbol)
                for symbol, _ in (
                    deserialize_replacements(
                        replacements
                    )
                )
            )

    for condition_solver in (
        solver.get(
            "conditions",
            {}
        ).values()
    ):
        result.update(
            _collect_cse_symbols(
                condition_solver
            )
        )

    return result


def assert_cse_temporaries_disjoint(
    solvers,
):
    """
    Mixed analytical/numerical solver blocks must not share temporary CSE
    variable names.
    """

    used = set()

    for solver in solvers:

        current = _collect_cse_symbols(
            solver
        )

        overlap = used.intersection(
            current
        )

        assert not overlap, (
            "CSE temporary symbols leaked across "
            f"solver blocks: {sorted(overlap)}"
        )

        used.update(current)


def assert_cse_dependency_order(
    solver,
):
    """
    Verify that each CSE replacement references only temporaries already
    defined earlier in its execution region.
    """

    for region_name, serialized_replacements in (
        solver.get("cse", {}).items()
    ):

        replacements = deserialize_replacements(
            serialized_replacements
        )

        defined = set()

        for temporary, expression in replacements:

            referenced_temporaries = {
                str(symbol)
                for symbol in expression.free_symbols
                if str(symbol).startswith(
                    "__ode_cse_"
                )
            }

            undefined = (
                referenced_temporaries
                - defined
            )

            assert not undefined, (
                f"CSE region {region_name} contains "
                "forward temporary dependencies: "
                f"{sorted(undefined)}"
            )

            defined.add(
                str(temporary)
            )
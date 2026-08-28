#
# test_common_subexpression_elimination.py
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

""""
This script provides a test of whether the original equation and CSE-optimised equation is preserved (==0). 
"""

import json 
import sympy 
from odetoolbox.expression_optimisation import (common_subexpression_elimination, restore_cse_expression, make_locals)

def test_cse_preserves_expression():

    x, y, a, b = sympy.symbols( # define sympy symbols 
        "x y a b",
        real=True
    )

    expressions = {
        "eq1": (x + y) * a,  # unoptimised eq
        "eq2": (x + y) * b,
    }

    replacements, reduced = (common_subexpression_elimination(expressions)) # perform cse 

    for name in expressions: # looping over per eq

        restored = restore_cse_expression(
            reduced[name],
            replacements
        )

        # simplify() is a SymPy function that rewrites mathematical expressions into a shorter or simpler form.
        difference = sympy.simplify(restored - expressions[name])

        # The text after the comma only shows up if difference != 0
        assert difference == 0, f"CSE failed. The restored expression differed from the original by: {difference}"

def test_cse_restores_operations_json():

    before = "/p/project1/paj2623/gray2/ode-toolbox-cse/diagnostics/outputs/amat_solver_baseline.json"
    after = "/p/project1/paj2623/gray2/ode-toolbox-cse/diagnostics/outputs/amat_solver_cse.json"

    with open(before, "r", encoding="utf-8") as file:
        model_before = json.load(file)[0]

    with open(after, "r", encoding="utf-8") as file:
        model_after = json.load(file)[0]

    locals_dict = make_locals(model_before)

    # Convert JSON CSE definitions into SymPy replacement tuples
    replacements = []

    for item in model_after["cse"]["propagators"]:
        temporary = sympy.Symbol(item["symbol"])
        replacement_expression = sympy.sympify(
            item["expression"],
            locals=locals_dict,
        )

        replacements.append(
            (temporary, replacement_expression)
        )

    # Now restore each propagator expression
    restored_propagators = {}

    for name, expression in model_after["propagators"].items():
        expression = sympy.sympify(
            expression,
            locals=locals_dict,
        )

        restored_propagators[name] = restore_cse_expression(
            expression,
            replacements,
        )

    # Compare with baseline
    for name, expression in model_before["propagators"].items():

        original = sympy.sympify(
            expression,
            locals=locals_dict,
        )

        restored = restored_propagators[name]

        difference = sympy.simplify(restored - original)

        assert difference == 0, (
            f"CSE failed for propagator {name}. "
            f"Difference: {difference}"
        )
        print("CSE Json file has been successfully restored. Mathematics has not been changed. ")
#
# test_cse_reduced_operation_count.py
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
This script provides tests that cse reduces the number mathematical operations from the original eq. 
"""

import sympy 
from odetoolbox.expression_optimisation import (common_subexpression_elimination, count_cse_operations, count_operations, parse_expression, make_locals)
import json 

def test_cse_preserves_expression():

    x, y, a, b = sympy.symbols( # define sympy symbols 
        "x y a b",
        real=True
    )

    expressions = {
        "eq1": (x + y) * a,  # unoptimised eq
        "eq2": (x + y) * b,
    }

    before = count_operations(expressions.values()) # count number of mathematical operations for original eq

    replacements, reduced = (common_subexpression_elimination(expressions)) # perform cse 

    after = count_cse_operations(replacements, reduced) # count number of mathematical operations for cse eq 

    # Print blocks for manual validation (visible pytest -s)
    print(f"\n Operations Before CSE: {before}")
    print(f" Operations After CSE: {after}")

    # enforce logic that optimisation must decrease operations 
    assert after < before, f"Optimisation failed, cost did not decrease. Before: {before}. After: {after}"


def test_cse_functional_json():

    before = ("/p/project1/paj2623/gray2/ode-toolbox-cse/diagnostics/outputs/amat_solver_baseline.json")
    after = ("/p/project1/paj2623/gray2/ode-toolbox-cse/diagnostics/outputs/amat_solver_cse.json")


    # Load baseline
    with open(before, "r", encoding="utf-8") as file:
        model_before = json.load(file)[0]


    # Load CSE result
    with open(after, "r", encoding="utf-8") as file:
        model_after = json.load(file)[0]

    # Create parser symbols
    locals_dict = make_locals(model_before)


    # Baseline cost
    before_propagator_cost = count_operations(
        model_before["propagators"].values(),
        locals_dict,
    )

    before_update_cost = count_operations(
        model_before["update_expressions"].values(),
        locals_dict,
    )

    before_total = (
        before_propagator_cost
        + before_update_cost
    )

    # CSE propagator cost
    after_propagator_cost = count_cse_operations(
        model_after["cse"]["propagators"],
        model_after["propagators"],
        locals_dict,
    )

    # CSE update-expression cost
    after_update_cost = count_cse_operations(
        model_after["cse"]["update_expressions"],
        model_after["update_expressions"],
        locals_dict,
    )

    after_total = (
        after_propagator_cost
        + after_update_cost
    )

    # --------------------------------------------------
    # Report
    # --------------------------------------------------

    print("\n=== Operation Count ===")

    print(
        f"Baseline propagators: "
        f"{before_propagator_cost}"
    )

    print(
        f"Baseline updates: "
        f"{before_update_cost}"
    )

    print(
        f"Baseline total: "
        f"{before_total}"
    )

    print()

    print(
        f"CSE propagators: "
        f"{after_propagator_cost}"
    )

    print(
        f"CSE updates: "
        f"{after_update_cost}"
    )

    print(
        f"CSE total: "
        f"{after_total}"
    )

    print()

    reduction = before_total - after_total
    percentage = (
        100 * reduction / before_total
    )

    print(
        f"Reduction: "
        f"{reduction}"
    )

    print(
        f"Reduction percentage: "
        f"{percentage:.2f}%"
    )

    # Actual test criterion
    assert after_total < before_total
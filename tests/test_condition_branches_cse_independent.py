#
# test_condition_branches_cse_independent.py
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

import sympy
import odetoolbox
from odetoolbox.expression_optimisation import (apply_cse_to_solver, restore_cse_expression)
import json 


"""
Independent Branch CSE: It ensures that when 'optimise_condition_branches=True', 
CSE is applied to both the default track AND special singularity branches.
Namespace Isolation: It proves that temporary CSE variables generated in one branch 
(e.g., __ode_cse_cond_0_prop_0) are completely disjoint from other branches, 
preventing variable leakage.
 Mathematical Correctness: It restores the optimized expressions back to their original 
algebraic forms using 'restore_cse_expression' and simplifies them to guarantee that 
the optimization pass did not mutate or alter the underlying physics equations.
"""


def test_condition_branches_cse_independently():

    x, y, a, b = sympy.symbols(
        "x y a b",
        real=True,
    )



    # Construct a mock solver profile containing a default region and a singularity region.
    # Note that both regions reuse the exact same subexpressions (x + y) and (x - y).
    solver = {
        "solver": "analytical",

        "conditions": {

            "default": {
                "update_expressions": {"x": a * (x + y), "y": b * (x + y)}},

            "(tau_1 == tau_2)": {
                "update_expressions": {"x": a * (x - y), "y": b * (x - y)}},
        }
    }

    # execute region cse optimisation pass 
    result = apply_cse_to_solver(solver, optimise_condition_branches=True)

    # Isolate the processed output regions
    default_branch = (result["conditions"]["default"]) # default solver with no known singularities 
    singularity_branch = (result["conditions"]["(tau_1 == tau_2)"]) # singularity branch 

    # test cse actually triggered and extracted terms 
    assert "cse" in default_branch
    assert "cse" in singularity_branch

    # [(temporary_var, original_expression)]
    default_replacements = (
        default_branch["cse"]
        ["update_expressions"]
    )

    singular_replacements = (
        singularity_branch["cse"]
        ["update_expressions"]
    )

    assert default_replacements
    assert singular_replacements

    # Collect only the left-hand names of the replacements (the new variable names)
    default_symbols = {temporary for temporary, _ in default_replacements}
    singular_symbols = {temporary for temporary, _ in singular_replacements}

    # No variable name generated in the default branch exists in the singularity branch.
    # If this fails, the downstream C++ generator will overwrite variables and cross wires.
    assert default_symbols.isdisjoint(singular_symbols)

    for condition, original_branch in (
        
        solver["conditions"].items()):

        optimized_branch = (result["conditions"][condition])

        replacements = (optimized_branch.get("cse", {}).get("update_expressions",[]))

        for variable, original in (
            original_branch["update_expressions"].items()):

            reduced = (optimized_branch["update_expressions"][variable])
            restored = restore_cse_expression(reduced,replacements,)

            assert sympy.simplify(restored - original) == 0


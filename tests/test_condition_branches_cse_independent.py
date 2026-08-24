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
from odetoolbox.expression_optimisation import (apply_cse_to_solver)
import json 

"""
specific test verifies non-finite expression check actually works, essentially as an crash protection test. 


"""


def test_condition_branches_cse_independently():

    x, y, a, b = sympy.symbols(
        "x y a b",
        real=True,
    )

    solver = {
        "solver": "analytical",

        "conditions": {

            "default": {
                "update_expressions": {"x": a * (x + y), "y": b * (x + y)}},

            "(tau_1 == tau_2)": {
                "update_expressions": {"x": a * (x - y), "y": b * (x - y)}},
        }
    }

    result = apply_cse_to_solver(solver, optimise_condition_branches=True)

    default_branch = (result["conditions"]["default"]) # default solver with no known singularities 
    singularity_branch = (result["conditions"]["(tau_1 == tau_2)"]) # singularity branch 

    assert "cse" in default_branch
    assert "cse" in singular_branch

    default_replacements = (
        default_branch["cse"]
        ["update_expressions"]
    )

    singular_replacements = (
        singular_branch["cse"]
        ["update_expressions"]
    )

    assert default_replacements
    assert singular_replacements
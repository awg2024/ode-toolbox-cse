#
# test_numeric_solver_cse.py
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
This script provides a test for numeric solver cse, can it optimise the RHS correctly,
update expressions, preserve its structure and mathematics? 
"""

import sympy 
from odetoolbox.expression_optimisation import apply_cse_to_solver, restore_cse_expression


def test_numeric_solver_cse():

    # define x,y as real sympy objects 
    x,y = sympy.symbols("x y", real=True)

    common_term = sympy.exp(x + y) # define common term as exp^x+y sympy object 

    # mocking solver data structure, in the expected format for sympy 
    solver_dict = {
        "solver": "numeric",
        "state_variables": ["x", "y"],
        "update_expressions": {
            "x": x + common_term,
            "y": y + 2 * common_term,
        },
    }

    original_expressions = dict(solver_dict["update_expressions"]) # original update expressions 

    result = apply_cse_to_solver(solver_dict) # apply function 

     
    assert result["solver"] == "numeric", "Failed: The solver type mutated or was lost."
    assert "cse" in result, "Failed: 'cse' key dictionary was never initialised."
    assert "update_expressions" in result["cse"], ("Failed: 'update_expressions' was not processed or missing inside the inner cse tracker.")

    replacements = result["cse"]["update_expressions"]

    assert len(replacements) > 0 # checking optimisation occured

    # check cse didnt remove eq
    assert set(result["update_expressions"]) == {"x", "y"}

    # prove equilance to original eq
    for variable, original in original_expressions.items():
        
        reduced = result["update_expressions"][variable] # Used singular 'variable'

        # Convert the reduced string expression back to a SymPy object for testing
        restored = restore_cse_expression(reduced, replacements) 

        assert sympy.simplify(restored - original) == 0


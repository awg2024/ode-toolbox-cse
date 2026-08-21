#
# test_analysis_numeric_cse.py
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
tests full pipeline of numeric cse analysis 
"""


from .context import odetoolbox
import sympy 
from odetoolbox.expression_optimisation import apply_cse_to_solver, restore_cse_expression


def test_analysis_numeric_cse():

    # mocking a .json from a model description file 
    indict = {
        "dynamics": [ # 
            {
                "expression":
                    "x' = x + exp(x + y)", # nonlinear 
                    "initial_value":
                        "0",
            },
            {
                "expression":
                    "y' = y + 2 * exp(x + y)", # nonlinear 
                    "initial_value":
                        "0",
            }, 
        ]    
    }

    result = odetoolbox.analysis( # calling the analysis primary function from init 
        indict,
        disable_analytic_solver=True, # forces entire calculation layout using a numerical integration path
        disable_stiffness_check=True, # disables GSL matrix operations to keep execution quick 
        enable_cse=True,  # switches on optimisation
    )

    # verification of dictionary and structure 
    assert len(result) == 1
    solver = result[0]
    assert solver["solver"] == "numeric"
    assert "update_expressions" in solver
    assert "cse" in solver
    assert "update_expressions" in solver["cse"]

    # verification of optimisation and serialisation 
    replacements = solver["cse"]["update_expressions"]
    assert len(replacements) > 0 # ensure that optimisation has occured 
    
    # To this (Index 0 Symbol name, Index 1 Math expression) - ensuring serialisation 
    assert isinstance(replacements[0][0], str), "Failed: The replacement symbol is not a string!"
    assert isinstance(replacements[0][1], str), "Failed: The replacement expression is not a string!"
    assert isinstance(solver["update_expressions"]["x"],str)
    assert isinstance(solver["update_expressions"]["y"],str)



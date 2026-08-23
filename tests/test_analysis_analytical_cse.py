
#
# test_analysis_analytical_cse.py
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
This script provides a test for an analytical solver cse, utilising the generate_propagator_solver()



"""

import sympy 
from odetoolbox.expression_optimisation import apply_cse_to_solver, restore_cse_expression
from .context import odetoolbox
def test_analysis_analytical_cse():

    indict = {
        "dynamics": [ # mocking indict that will be passed into _analysis 
            { 
                "expression": "x' = -x / tau",
                "initial_value": "1",
            },
            {
                "expression": "y' = -y / tau",
                "initial_value": "2",
            }
        ]
    }

    result = odetoolbox.analysis(
        indict, 
        disable_stiffness_check=True,
        disable_singularity_detection=True, 
        enable_cse=True
    )

    # verification of dictionary and structure 
    assert len(result) == 1
    solver = result[0]
    assert solver["solver"] == "analytical"
    assert "propagators" in solver
    assert "cse" in solver
    assert "propagators" in solver["cse"]

    # verification of optimisation and serialisation 
    replacements = solver["cse"]["propagators"]
    assert len(replacements) > 0 # ensure that optimisation has occured 
    
    # To this (Index 0 Symbol name, Index 1 Math expression) - ensuring serialisation 
    replacement = replacements[0]

    assert isinstance(replacement["symbol"], str)

    assert isinstance(replacement["expression"], str)

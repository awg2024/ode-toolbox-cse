
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

def test_analytical_solver_cse():

    indict = {
        "dynamics": [
            {
                "expression":
                    "x' = -x / tau",
                "initial_value":
                    "1",
            },
            {
                "expression":
                    "y' = -y / tau",
                "initial_value":
                    "2", 
            }
        ]
    }

    result = odetoolbox.analysis(
        indict, 
        disable_stiffness_check=True,
        disable_singularity_detection=True, # disable because analytical produces a conditional structure in the json 
        enable_cse=True)

    assert len(result) == 1 

    solver = result[0]

    assert solver["solver"] == "analytical"
    assert "propagators" in solver
    assert "cse" in solver
    assert "propagators" in solver["cse"]

    propagator_replacements = (solver["cse"]["propagators"])

    assert len(propagator_replacements) > 0 

    first_replacement = (propagator_replacements[0])

    # stored as replacement symbol and expression tuple access 
    assert isinstance(first_replacement[0], str), "Failed: The replacement symbol name is not a string!"
    assert isinstance(first_replacement[1], str), "Failed: The replacement expression body is not a string!"


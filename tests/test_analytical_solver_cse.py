#
# test_analytical_solver_cse.py
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


def test_analytical_solver_cse():

    
    # This defines the standard physical time step parameter (__h) and an exponential decay time constant (tau).
    h, tau = sympy.symbols("__h tau", real=True) 
    common_propagator = sympy.exp(-h / tau) # common propagator found across neuronal models 

    P_x = sympy.Symbol("__P__x__x", real=True) # defining propagators 
    P_y = sympy.Symbol("__P__y__y", real=True)

    solver_dict = { # mocking an analytical matrix layout 
        "solver": "analytical",
        "state_variables": [
            "x",
            "y",
        ],

        "propagators": {
            "__P__x__x": common_propagator, # defining the RHS as the common propagator 
            "__P__y__y": common_propagator,
        },

        "update_expressions": {
            "x": P_x * sympy.Symbol("x", real=True), # x,y real sympy objects
            "y": P_y * sympy.Symbol("y", real=True),
        },
    }

    original_propagators = dict(solver_dict["propagators"]) # saving benchmark 

    result = apply_cse_to_solver(solver_dict) # apply function 

     
    assert result["solver"] == "analytical", "Failed: The solver type mutated or was lost."
    
    assert "propagators" in result["cse"]

    replacements = result["cse"]["propagators"]

    assert len(replacements) > 0 # checking optimisation occured

    # prove equilance to original eq
    for variable, original in original_propagators.items(): # very similar code to numeric test 
        
        reduced = result["propagators"][variable] # Used singular 'variable'

        # Convert the reduced string expression back to a SymPy object for testing
        restored = restore_cse_expression(reduced, replacements) 

        assert sympy.simplify(restored - original) == 0


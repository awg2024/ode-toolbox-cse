#
# cse_utils.py
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

import pytest
import sympy 
import odetoolbox.expression_optimisation as eo 
from .context import odetoolbox

"""
A related script cse_utils.py can be found which contains functions for the testing of CSE. 
"""

def _run_cse_analysis(indict, *, force_numeric=False):
    
    kwargs = { # define keyword arguments 
        "disable_singularity_detection": True,
        "disable_stiffness_check": True,
        "log_level": "DEBUG" 
    }

    if force_numeric:
        kwargs["disable_analytic_solver"] = True 

    baseline, _, _ = odetoolbox._analysis(indict, enable_cse=False, **kwargs) # disabled 

    optimised, shape_sys, shapes = odetoolbox._analysis(indict, enable_cse=True, **kwargs) # enabled 



def restore_cse_expression(expression, replacements):
    # reversal of cse back to the original form.

    restored = expression # sypy expressions are immutable, each subs() returns a transformed expression

    for temporary, replacement_expression in reversed(replacements): # reversed as later temp variables depends on earlier definitions
        restored = restored.subs(temporary, replacement_expression) #.subs finds temporary variable and swaps it for it's actual mathematical replacement

    return restored

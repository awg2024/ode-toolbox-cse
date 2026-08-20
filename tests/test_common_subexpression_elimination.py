#
# XXXX.py
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
This script provides a simple subexpression eliminiation (early dev)

After CSE implementation we need tests for the following: 
    • analytical solver CSE;
    
    • numerical solver CSE;
    
    • conditional singularity solvers;
    
    • expression serialization;
    
    • nested CSE temporaries;
    
    • CSE disabled;
    
    • unchanged legacy JSON output when disabled.
"""

import sympy 
from odetoolbox.expression_optimisation import (common_subexpression_elimination)

def test_basic_common_subexpression_elimination():

    x, y, a, b = sympy.symbols( # define sympy symbols 
        "x y a b",
        real=True
    )

    expressions = {
        "eq1": (x + y) * a,  # unoptimised eq
        "eq2": (x + y) * b,
    }

    replacements, reduced = (common_subexpression_elimination(expressions))

    assert len(replacements) == 1

    temporary, temporary_expression = replacements[0]

    # checks, we know the expected cse outcome so we can compare. 
    assert temporary_expression == x + y
    assert reduced["eq1"] == temporary * a
    assert reduced["eq2"] == temporary * b
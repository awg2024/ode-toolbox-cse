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
from odetoolbox.expression_optimisation import (common_subexpression_elimination, count_cse_operations, count_operations)

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



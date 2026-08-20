#
# expression_optimisation.py
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
system_of_shapes.py - mathematical solver
sympy_helpers.py - currently contains low-level utilities

This script provides a transformational pass for sympy expressions evaluating: 
- cost calculations
- dependency analysis
- temporary validation
- operation counting
- future alegbraic optimisations...

"""



import sympy 

def common_subexpression_elimination(expressions, symbol_prefix="__ode_cse_tmp__"):
    """
    custom wrapper to perform common subexpression elimination across mapping of a named sympy expression
    """

    if not expressions: # if expressions are empty 
        return [], {}
    
    #match values with names
    expressions_names = list(expressions.keys()) 
    expressions_values = [expressions[name] for name in expressions_names] # sympy accepts ordered list of mathematical expressions 
    # splitting variable names and raw sympy math trees 

    temporary_symbols = sympy.numbered_symbols(symbol_prefix) # infinite generator for the temp local variables to hold isolated math subexpressions 

    # perform cse
    replacements, reduced_values = sympy.cse(
        expressions_values, # 
        symbols=temporary_symbols,
        optimizations=None, 
        order="canonical" # forces ordering deterministically 
    )

    reduced_expressions = dict(zip(expressions_names, reduced_values))

    return replacements, reduced_expressions

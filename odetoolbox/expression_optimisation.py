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

def restore_cse_expression(expression, replacements):
    # reversal of cse back to the original form. 

    restored = expression # creates copy 

    for temporary, replacements in reversed(replacements): # reversed as later temp depends on earlier 
        restored = restored.subs(temporary, replacements) #.subs finds temporary variable and swaps it for it's actual mathematical replacement 
    
    return restored 


def count_operations(expressions):

    return sum( # counting number of mathematical operations from the original equation 
        int(sympy.count_ops(expression))
        for expression in expressions
    )

def count_cse_operations(replacements, reduced_expressions): 

    # count operations in the main simplified equations
    reduced_cost = sum(int(sympy.count_ops(expr)) for expr in reduced_expressions.values())

    # Count operations inside the temporary placeholder variables
    replacement_cost = sum(int(sympy.count_ops(expr)) for _, expr in replacements)

    return replacement_cost + reduced_cost # total final cost of expressions and temporary values 


def apply_cse_to_solver(solver, symbol_prefix="__ode_cse_"): # prefix to use downstream of this value belongs to a cse optimisation, enabling a clean variable handling
    """
    It's import to understand ODE toolbox output structure (generate_solver_dict_based_on_propagator_matrix(), generate_numeric_solver())
    If we extract out temporary values from the entirity of nestml at once we would lead to crashes in the ODEs since these equations may be updated differently 
    and at different steps therefore the temporary value would have no where to live. 
    we will have different regions of cse extraction and optimisation: update expression, propagator expression, singularity-condition branch. 
    """

    # if "propagators" in solver etc.. should we have a flag? 

    if "update_expressions" in solver:
        
        replacements, reduced = (common_subexpression_elimination(solver["update_expressions"], symbol_prefix="__ode_cse_update"))

        result["update_expressions"] = reduced

    # for a singularity-expression, they will have their own propagators and update_expressions

    return result 


def serialize_replacements(replacements):
"""
serializing the SymPy varible into a raw text string inside a json file so that code generator can read it and write it out as a line of C++ code 
"""
    return[ 
        {
            "symbol": str(symbol),
            "expression": str(expression),
        }
        for symbol, expression in replacements 
    ]


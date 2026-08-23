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
    custom wrapper to perform common subexpression elimination across mapping of a 
    named sympy expression
    """

    if not expressions: # if expressions are empty 
        return [], {}
    
    # expression names in specific order 
    expressions_names = list(expressions.keys()) 

    # sympy accepts ordered list of mathematical expressions in same order as expression_names
    expressions_values = [expressions[name] for name in expressions_names] 

    # check for sympy objects 
    if not all(isinstance(expression, sympy.Basic) for expression in expressions_values):
        raise TypeError("CSE expects Sympy objects. String serialisation has not occured. ")

    #  infinite generator for the temp local variables to hold isolated math subexpressions 
    temporary_symbols = sympy.numbered_symbols(symbol_prefix) 

    # perform cse
    replacements, reduced_values = sympy.cse( 
        expressions_values, # 
        symbols=temporary_symbols,
        optimizations=None, 
        order="canonical" # forces ordering deterministically (e.g., A relies on B)
    )

    # fuse cse expressions with values, maintain order of eq. 
    reduced_expressions = dict(zip(expressions_names, reduced_values))

    # replacement is returned as tuples where each tupe contains the temp var_name with the maths block
    return replacements, reduced_expressions

def restore_cse_expression(expression, replacements):
    # reversal of cse back to the original form. 

    restored = expression # creates copy 

    for temporary, replacement_expression in reversed(replacements): # reversed as later temp depends on earlier 
        restored = restored.subs(temporary, replacement_expression) #.subs finds temporary variable and swaps it for it's actual mathematical replacement 
    
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
    Apply CSE independently to the different execution regions of an ODE-toolbox solver (e.g., update, propagators)

    All expressions remain SymPy objects here. Serialization happens later in _analysis() for a safe JSON-safe solver dict. 
    """

    result = dict(solver)

    cse_data = {}
    
    # analytical propagator solver region (non-linear ODE)
    if "propagators" in solver: 
        
        replacements, reduced = common_subexpression_elimination(
            solver["propagators"],
            symbol_prefix=symbol_prefix + "prop_")

         # keep sympy object here 
        result["propagators"] = reduced

        if replacements:
            cse_data["propagators"] = (replacements)

    # numerical state update region (linear ODE)  
    if "update_expressions" in solver:

        replacements, reduced = common_subexpression_elimination(
            solver["update_expressions"],
            symbol_prefix=symbol_prefix + "update_")
        
        # keep sympy object here 
        result["update_expressions"] = reduced

        if replacements:
            cse_data["update_expressions"] = (replacements)

    # Only add CSE metadata when something was actually extracted.
    if cse_data:
        result["cse"] = cse_data
    
    
    return result



def serialize_replacements(replacements):
    """
    Convert CSE replacement tuples to JSON-safe data
    """

    return [{
            "symbol": str(symbol),
            "expression": str(expression)}
            for symbol, expression in replacements]


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
Expression optimisation script for the common subexpression elimination (CSE) flag 
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


def apply_cse_to_solver(solver, symbol_prefix="__ode_cse_", optimise_condition_branches=False): 
# prefix to use downstream of this value belongs to a cse optimisation inside json

    """
    Apply CSE independently to the different execution regions of an ODE-toolbox solver (e.g., update, propagators, singularity-conditions)
    Nested replacements serialization is performed safely at the boundaries of each region.
    """

    result = dict(solver)
    
    cse_data = {}
    
    if "conditions" in solver:

        # enable conditions to pass through depending on flag
        if not optimise_condition_branches:
            return dict(solver)
        
        result = dict(solver)
        optimised_conditions = {}

        # wrap condition, branch so python understand how to unpack a sub-tuple 
        for branch_index, (condition, branch) in enumerate(solver["conditions"].items()):

            branch_prefix = (symbol_prefix + f"cond_{branch_index}") # distinct condtion blocks get their own unique prefix

            # count mathematical cost before cse 
            raw_expressions = [] # count number of operations
            if "update_expressions" in branch:
                raw_expressions.extend(branch["update_expressions"].values())
            if "propagators" in branch:
                raw_expressions.extend(list(branch["propagators"]))
            before_cost = count_operations(raw_expressions)

            # apply localised region cse optimisation 
            optimized_region = _apply_cse_to_expression_region(branch, symbol_prefix=branch_prefix)

            # Loop through your active target regions 
            for region_key in ["update_expressions", "propagators"]:
                branch_replacements = optimized_region["cse"][region_key]
                branch_reduced = optimized_region[region_key]

                # Calculate your true optimized operational weight
                after_cost = count_cse_operations(branch_replacements, branch_reduced)
            
            # if the branch optimisation didn't reduced cost 
            if (replacements and after_cost < before_cost):
                optimised_conditions[condition] = dict(branch)
                continue  # Skip serialization and proceed to the next condition branch
            
            if "cse" in optimized_region:
                optimized_region["cse"] = _serialize_replacements_metadata(optimized_region["cse"])   # Serialize the metadata after your cost metrics have been calculated

        # Store the completely optimized and validated block back into your conditions matrix
        optimised_conditions[condition] = optimized_region


    #  analytical propagator solver region (ordinary solver)
    if "propagators" in solver:  # searching high-level dict 
        
        # flatten matrix and calculate original propagator mathematical cost
        before_cost = count_operations(list(solver["propagators"]))

        replacements, reduced = common_subexpression_elimination(
            solver["propagators"],
            symbol_prefix=symbol_prefix + "prop_")

        after_cost = count_cse_operations(replacements, reduced)

        if replacements and after_cost >= before_cost:
            result["propagators"] = solver["propagators"]  # Revert reduced expression back to the untouched original matrix layout

         # keep sympy object here 
        result["propagators"] = reduced

        if replacements:
            cse_data["propagators"] = (replacements)

    # ------- numerical state update region (ordinary solver)  
    if "update_expressions" in solver: # searching high level dict 


        before_cost = count_operations(list(solver["update_expressions"]))

        replacements, reduced = common_subexpression_elimination(
            solver["update_expressions"],
            symbol_prefix=symbol_prefix + "update_")

        after_cost = count_cse_operations(replacements, reduced)

        if replacements and after_cost >= before_cost:
            result["propagators"] = solver["propagators"]  # Revert reduced expression back to the untouched original matrix layout

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
    Convert CSE replacement tuples to JSON-safe data (high-level dicts: update_expressions, propagators)
    """

    return [{
            "symbol": str(symbol),
            "expression": str(expression)}
            for symbol, expression in replacements]


def _serialize_replacements_metadata(region):
    """
    serialize cse replacement tuples nested in a solver (e.g., when solving a condition json we need to search deeper than top-level dicts)
    """

    # pass if cse wasn't conducted
    if not "cse" in region:
        return
    
    for expression_region, replacements in (region["cse"].items()):

        # call the indivudal serialise functon once within the solver 
        region["cse"[expression_region] = serialize_replacements(replacements)]




def _apply_cse_to_expression_region(region, symbol_prefix):

    """
    Apply CSE inside one execution region of the system_of_shapes.generate_propagator_solver(). 
    The condition controlling execution region is not modified, and therefore this function should never 
    receive multiple singularity branches as this function does not hold logic for interpreting these singularities. 
    """
    
    result = dict(region)
    cse_data = {}  # use a dict instead of a list [] to avoid KeyError downstream

    # collect all expressions safely into a list for the non-finite check
    all_math_expressions = []

    if "propagators" in region and region["propagators"] is not None:   # searching for propagator expression inside this conditional sympy object dict
        # If it's a SymPy Matrix or list, extend your tracking bucket safely
        all_math_expressions.extend(list(region["propagators"]))
        
    if "update_expressions" in region and region["update_expressions"] is not None:     # searching for update expressions inside this conditional sympy object dict 
        all_math_expressions.extend(list(region["update_expressions"]))
    

    # safety check before running cse 
    if _contains_nonfinite_expression(all_math_expressions): # second sanity check for singularities
        print(f"Non-finite/Infinity expression detected inside region prefix: {symbol_prefix}")
        return result

    # safety check before running cse 
    if _contains_internal_control_flow(all_math_expressions):
        print(f"Hidden nest Sympy.Piecewise logic {symbol_prefix}")
        return result

    if "propagators" in region and region["propagators"] is not None:
        replacements, reduced = common_subexpression_elimination( # perform cse on propagators inside condition
            region["propagators"], 
            symbol_prefix=(symbol_prefix + "prop_")
        )
        if replacements: 
            result["propagators"] = reduced
            cse_data["propagators"] = replacements 
    
    if "update_expressions" in region and region["update_expressions"] is not None: # perform cse on update_expressions inside condition
        replacements, reduced = common_subexpression_elimination(
            region["update_expressions"], 
            symbol_prefix=(symbol_prefix + "update_")
        )
        if replacements: 
            result["update_expressions"] = reduced
            cse_data["update_expressions"] = replacements 
    
    if cse_data:
        result["cse"] = cse_data # output data 

    return result


    # searching for propagator expression inside this conditional sympy object dict
    if "propagators" in region:

        replacements, reduced = common_subexpression_elimination(region["propagators"], symbol_prefix=(symbol_prefix + "prop_"))

        if replacements: 
            result["propagators"] = reduced
            cse_data["propagators"] = replacements 
    
    # searching for update expressions inside this conditional sympy object dict 
    if "update_expressions" in region:

        replacements, reduced = common_subexpression_elimination(region["update_expressions"], symbol_prefix=(symbol_prefix + "update_"))

        if replacements: 
            result["update_expressions"] = reduced
            cse_data["update_expressions"] = replacements 

    if cse_data:
        result["cse"] = cse_data

    return result 


def _contains_nonfinite_expression(expressions):

    """
    This is check before CSE occurs to check for symbolic infinities. This looks at an equation and determiens
    that it will always evaluate to infinity or divide by 0, regardless of the numerical values you pass. 
    This is secondary sanity check, as all symbolic infinities should theortically be filtered out 
    by the singularity conditions. 
    """

    invalid_values = (
    sympy.zoo, # complex infinity 
    sympy.oo, # infinity
    -sympy.oo, # negative infinity 
    sympy.nan # NaN
    )

    for expression in expressions:
        if hasattr(expression, "has"):
            for invalid in invalid_values:
                if expression.has(invalid):
                        return True # if the value contains these invalid values return true 
                    
    return False


def _contains_internal_control_flow(expressions):

    """
    If a sympy.Piecewise object is hidden inside an expression, it introduces hidden branching logic. If CSE blindly pulls an equation out from inside a 
    Piecewise condition and places it at the global scope, it forces the CPU to compute it all the time. This ruins your conditional optimization and can 
    lead to runtime NaN crashes or division-by-zero errors. This acts a secondary safety check before cse. 
    """

    # If passed a dictionary (like your solver_dict["update_expressions"])
    if isinstance(expressions, dict):
        return any(expression.has(sympy.Piecewise) for expression in expressions.values())
    
    
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
Expression optimisation helpers used by the ODE-toolbox analysis pipeline. 

This module contains only runtime common subexpression elimination (CSE) functionality. Test-specific 
reconstruction and validation helpers live in tests/cse_test_utils.py 
"""

import logging 
import sympy 

def common_subexpression_elimination(expressions, symbol_prefix="__ode_cse_tmp__"):
    """
    custom wrapper to perform common subexpression elimination across mapping of a
    named sympy expression while preserving expression ordering. 
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


def count_operations(expressions):
    """
    Count the total symbolic operations in an iterable of SymPy expressions. 
    """
    
    #counting number of mathematical operations from the original equation 
    return sum(int(sympy.count_ops(expression)) for expression in expressions)


def count_cse_operations(replacements, reduced_expressions): 
    """
    Count operations required after CSE. helper includes calculating the extracted temporary expressions & reduced output expressions 
    """

    # Extract operational cost for main simplified expressions
    reduced_cost = sum(int(sympy.count_ops(expr)) for expr in reduced_expressions.values())

    # count operations inside the temporary placeholder variables
    replacement_cost = sum(int(sympy.count_ops(expr)) for _, expr in replacements)

    # assign an assignment/lookup penalty for every subexpression extracted
    temporary_overhead = len(replacements)
    
    return (replacement_cost + reduced_cost + temporary_overhead) # total final cost of expressions and temporary values and conservative overhead memory 



def _apply_cse_to_solver(solver, symbol_prefix="__ode_cse_", optimise_condition_branches=False):
    """
    Apply CSE independently to the different execution regions and retain it only when the symbolic operation count is reduced
    """

    if not expressions: 
        return [], expressions
    
    replacements, reduced = common_subexpression_elimination(expressions, symbolic_prefix=symbol_prefix)

    if not replacements:
        return [], expressions
    
    before_cost = count_operations(expressions.values())
    after_cost = count_cse_operations(replacements, reduced)
    
    if after_cost >= before_cost:
        return [], expressions # return original form 

    return replacements, reduced 



def _run_profitable_cse(expressions, symbol_prefix, solver_name="unknown", region_name="unknown"):

    """
    Apply CSE to one execution region and retain it only when the mathematical count of oeprations is reduced
    """

    if not expressions:
        return [], expressions
    
    replacements, reduced = (common_subexpression_elimination(expressions, symbol_prefix=symbol_prefix))

    logger = logging.getLogger(__name__)

    if not replacements:
        logger.debug("CSE [%s]: no common subexpression found. Replacements: [%s]", symbol_prefix, replacements)
        return [], expressions

    # mixed updated is returning [] !!! 
    
    original_ops = count_operations(expressions.values())
    replacement_ops = sum(int(sympy.count_ops(expr)) for _, expr in replacements)
    reduced_ops = sum(int(sympy.count_ops(expr)) for expr in reduced.values()) #  i want to use my count_cse_operations here? 
    temporary_count = len(replacements)
    temporary_penalty = temporary_count

    before_cost = original_ops
    after_cost = replacement_ops + reduced_ops + temporary_penalty

    if before_cost > 0:
        saving = before_cost - after_cost
        reduction_percentage = (saving / before_cost) * 100
    else:
        saving = 0
        reduction_percentage = 0.0

    decision = "ACCEPT" if after_cost < before_cost else "REJECT"

    logger.debug(
        "[CSE] solver=%s\n"
        "[CSE] region=%s\n\n"
        "original_expression_ops = %d\n"
        "replacement_ops         = %d\n"
        "reduced_expression_ops  = %d\n"
        "temporary_count         = %d\n"
        "temporary_penalty       = %d\n"
        "estimated_before         = %d\n"
        "estimated_after          = %d\n"
        "estimated_saving         = %d\n"
        "estimated_reduction      = %.2f%%\n"
        "decision                 = %s\n",
        solver_name, region_name,
        original_ops, replacement_ops, reduced_ops, temporary_count, temporary_penalty,
        before_cost, after_cost, saving, reduction_percentage, decision)

    if decision == "REJECT":
        logging.getLogger(__name__).debug(
            "CSE [%s] optimisation is rejected because it doesn't reduce symbolic operation count",
            symbol_prefix)
        return [], expressions

    logging.getLogger(__name__).debug(
        "CSE [%s]: accepted",
        symbol_prefix)
    
    return replacements, reduced

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
    Piecewise condition and places it at the global scope (when it has a local conditional specifications), it forces the CPU to compute it all the time. 
    
    This ruins your conditional optimization and can lead to runtime NaN crashes or division-by-zero errors. This acts a secondary safety check before cse.
    """

    # returns True/False for any sympy Piecewise is present
    return any(isinstance(expression, sympy.Basic) and expression.has(sympy.Piecewise) for expression in expressions)



def _apply_cse_to_expression_region(region, symbol_prefix, solver_name="unknown"):

    """
    Apply CSE independently inside one execution region.

    The condition controlling execution region is not modified, and therefore this function should never
    receive multiple singularity branches as this function does not hold logic for interpreting these singularities.

    'propagator' 'update_expressions' are also handled differently as they are executed in different contexts down stream. 
    """

    result = dict(region)
    cse_data = {}  # use a dict instead of a list [] to avoid KeyError downstream

    # collect all expressions safely into a list for the non-finite check
    all_math_expressions = []

    if region.get("propagators"): # collect update_expressions values inside all_math_expressions
        all_math_expressions.extend(region["propagators"].values())

    if region.get("update_expressions"): # collect propagator values inside all_math_expressions
        all_math_expressions.extend(region["update_expressions"].values())

    # safety check before running cse in singularity
    if _contains_nonfinite_expression(all_math_expressions): # second sanity check for singularities
        logger = logging.getLogger(__name__)
        logger.debug("Skipping CSE for region %s: non-finite symbolic expression detected", symbol_prefix)
        return result

    # safety check before running cse in singularity
    if _contains_internal_control_flow(all_math_expressions):
        logger = logging.getLogger(__name__)
        logger.debug("Skipping CSE for region %s: nested SymPy Piecewise expression detected", symbol_prefix)
        return result

    # analytical 
    if region.get("propagators"):
        replacements, reduced = _run_profitable_cse(
            region["propagators"], symbol_prefix + "prop_",
            solver_name=solver_name, region_name="propagators")
        if replacements:
            result["propagators"] = reduced
            cse_data["propagators"] = replacements

    # numerical state update expressions 
    if region.get("update_expressions"):
        replacements, reduced = _run_profitable_cse(
            region["update_expressions"], symbol_prefix + "update_",
            solver_name=solver_name, region_name="update_expressions")
        if replacements:
            result["update_expressions"] = reduced
            cse_data["update_expressions"] = replacements
            
    if cse_data:
        result["cse"] = cse_data # output data

    return result


def apply_cse_to_solver(solver, symbol_prefix="__ode_cse_", optimise_condition_branches=False):

    """
    Apply CSE to one ODE-toolbox solver dictionary. Singularity branches are treated as independent execution regions. 
    """

    result = dict(solver)
    solver_name = solver.get("solver", "unknown")

    if "conditions" in solver: # handle singularity conditions

        # enable conditions to pass through depending on flag
        if not optimise_condition_branches:
            return dict(solver)

        result = dict(solver)
        optimised_conditions = {}

        # wrap condition, branch so python understand how to unpack a sub-tuple
        for branch_index, (condition, branch) in enumerate(solver["conditions"].items()):

            branch_prefix = (symbol_prefix + f"cond_{branch_index}_") # distinct condition blocks get their own unique prefix

            optimised_conditions[condition] = _apply_cse_to_expression_region(branch, symbol_prefix=branch_prefix, solver_name=solver_name) # apply cse to each independent branch

        result["conditions"] = optimised_conditions

        return result # return result if we've conducted cse singularity

    # ordinary analytical or numerical solver passing. 
    return _apply_cse_to_expression_region(solver, symbol_prefix=symbol_prefix, solver_name=solver_name)


def _apply_cse_to_solver_blocks(solver_blocks, optimise_condition_branches=False):

    """
    Apply CSE independently to every solver block produced by ODEtoolbox 
    
    A mixed system can contain both analytical and numerical blocks. Give each block a seperate tmp variable namespace so that CSE
    temporaries cannot leak or collide across solver blocks. 
    """

    if not isinstance(solver_blocks, list):
        raise TypeError("_apply_cse_to_solver_blocks() expects a list of solver dictionaries",
                        f"received: {type(solver_blocks).__name__}")


    result = []
    multiple_blocks = len(solver_blocks) > 1 

    for block_index, solver in enumerate(solver_blocks):

        logging.getLogger(__name__).debug(
            "Applying CSE to solver block %d (%s)",
            block_index,
            solver.get("solver", "unknown"))

        if multiple_blocks:
            symbol_prefix = f"__ode_cse_solver_{block_index}_"
        else:
            symbol_prefix = f"__ode_cse_"
        
        result.append(apply_cse_to_solver(solver, symbol_prefix=symbol_prefix, optimise_condition_branches=(optimise_condition_branches)))

    return result 

def serialize_replacements(replacements):
    """
    Convert CSE replacement tuples into JSON-safe metadata. 
    """
    # import pdb;pdb.set_trace()

    result = {}
    for symbol, expr in replacements:
        result[str(symbol)] = str(expr)

    return result


def _serialize_replacements_metadata(region):
    """
    Serialize CSE replacement metadata belonging to one solver region.
    """

    # pass if cse wasn't conducted
    if not "cse" in region:
        return

    for expression_region, replacements in (list(region["cse"].items())):

        # call the indivudal serialise functon once within the solver
        region["cse"][expression_region] = serialize_replacements(replacements)

def _find_non_json_serializable(obj, path="root"):
    # Clear terminal conditions for valid JSON types
    if isinstance(obj, (str, int, float, bool, type(None))):
        return []

    # Handle dictionaries safely
    if isinstance(obj, dict):
        problems = []
        for key, value in obj.items():
            if not isinstance(key, (str, int, float, bool, type(None))):
                # Wrapped in a proper 3-element tuple inside the append
                problems.append((f"{path}.<key>", type(key).__name__, repr(key)))
            problems.extend(_find_non_json_serializable(value, path=f"{path}[{key!r}]"))
        return problems 
    
    # Handle sequences
    if isinstance(obj, (list, tuple)):
        problems = []
        for index, value in enumerate(obj):
            problems.extend(_find_non_json_serializable(value, path=f"{path}[{index}]"))
        return problems 
        
    # Catch-all for non-serializable objects (like SymPy Symbols)
    # Corrected to a uniform list containing a single 3-element tuple
    return [(path, type(obj).__name__, repr(obj))]

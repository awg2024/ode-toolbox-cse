import itertools
import time
import sympy
from amat_diagnostic_utils import (build_amat_system)
from odetoolbox.singularity_detection import (SingularityDetection)
from odetoolbox.sympy_helpers import (SymmetricEq)


"""
This script is utilised after the contigious and non-contigious checks investigating 
why would we get an explosion of 64 combinations singularity checks globally on a 8x8 matrix? 
Is there a problem with the generate_propagator_matrix function? or something upstream? 
"""


# helper function 
(model,shape_sys,analytic_sys,classification) = build_amat_system()

print("Generating one propagator")

start = time.perf_counter()

# generate propaga
P = (analytic_sys._generate_propagator_matrix(analytic_sys.A_))


print(
    "Propagator time:",
    time.perf_counter() - start,
)


print(
    "\nDetecting singularities..."
)


start = time.perf_counter()


conditions = (
    SingularityDetection
    .find_propagator_singularities(
        P,
        analytic_sys.A_,
    )
)


conditions = conditions.union(
    SingularityDetection
    .find_inhomogeneous_singularities(
        analytic_sys.A_,
        analytic_sys.b_,
    )
)


conditions = (
    SingularityDetection
    ._remove_duplicate_conditions(
        conditions
    )
)


conditions = list(
    conditions
)


print(
    "Detection time:",
    time.perf_counter() - start,
)


print(
    "\n======================================"
)

print(
    "DETECTED CONDITIONS"
)

print(
    "======================================"
)


for index, condition in enumerate(
    conditions
):

    print(
        index,
        ":",
        condition,
    )


N = len(
    conditions
)


print(
    "\nNumber of independent Boolean "
    "flags currently assumed:",
    N,
)


print(
    "Raw permutations:",
    2 ** N,
)


# ==================================================
# CONSISTENCY CHECK
# ==================================================

def condition_set_is_consistent(
    condition_set,
):

    parent = {}

    def find(item):

        if item not in parent:
            parent[item] = item

        if parent[item] != item:
            parent[item] = find(
                parent[item]
            )

        return parent[item]

    def union(left, right):

        root_left = find(
            left
        )

        root_right = find(
            right
        )

        if root_left != root_right:
            parent[
                root_right
            ] = root_left

    inequalities = []

    # First establish equality classes.
    for condition in condition_set:

        if isinstance(
            condition,
            SymmetricEq,
        ):

            union(
                condition.lhs,
                condition.rhs,
            )

        else:

            inequalities.append(
                (
                    condition.lhs,
                    condition.rhs,
                )
            )

    # Then check whether any inequality
    # contradicts those equalities.
    for left, right in inequalities:

        if find(left) == find(right):
            return False

    return True


consistent = 0
inconsistent = 0

consistent_nondefault = 0


for permutation in itertools.product(
    [False, True],
    repeat=N,
):

    condition_set = []

    for condition, holds in zip(
        conditions,
        permutation,
    ):

        if holds:

            condition_set.append(
                condition
            )

        else:

            condition_set.append(
                sympy.Ne(
                    condition.lhs,
                    condition.rhs,
                )
            )

    if condition_set_is_consistent(
        condition_set
    ):

        consistent += 1

        if any(
            isinstance(
                condition,
                SymmetricEq,
            )
            for condition
            in condition_set
        ):

            consistent_nondefault += 1

    else:

        inconsistent += 1


print(
    "\n======================================"
)

print(
    "BOOLEAN REGIME ANALYSIS"
)

print(
    "======================================"
)


print(
    "Raw regimes:",
    2 ** N,
)


print(
    "Consistent regimes:",
    consistent,
)


print(
    "Inconsistent regimes:",
    inconsistent,
)


print(
    "Consistent non-default "
    "alternate regimes:",
    consistent_nondefault,
)


# --------------------------------------------------
# Diagnostic only:
# What if __h == 0 is excluded?
# DO NOT use this to change production logic yet.
# --------------------------------------------------

conditions_without_h = [
    condition
    for condition in conditions
    if "__h"
    not in {
        str(symbol)
        for symbol
        in condition.free_symbols
    }
]


print(
    "\nConditions excluding __h:",
    len(
        conditions_without_h
    ),
)


print(
    "Raw regimes without __h:",
    2 ** len(
        conditions_without_h
    ),
)
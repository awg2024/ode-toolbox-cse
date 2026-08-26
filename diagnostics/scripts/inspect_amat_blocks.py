import numpy as np
import scipy.sparse

from amat_diagnostic_utils import (
    build_amat_system,
)

from odetoolbox.system_of_shapes import (
    get_block_diagonal_blocks,
    GetBlockDiagonalException,
)


(
    model,
    shape_sys,
    analytic_sys,
    classification,
) = build_amat_system()


# --------------------------------------------------
# Work specifically with the analytical subsystem A.
# --------------------------------------------------

A = np.array(
    analytic_sys.A_,
    dtype=object,
)


print(
    "\n======================================"
)

print(
    "ANALYTICAL SYSTEM MATRIX"
)

print(
    "======================================"
)

print(
    "A shape:",
    A.shape,
)


# --------------------------------------------------
# Construct the exact undirected connectivity graph
# used by get_block_diagonal_blocks().
# --------------------------------------------------

connectivity = (
    (A != 0)
    |
    (A.T != 0)
)


n_components, labels = (
    scipy.sparse.csgraph.connected_components(
        connectivity
    )
)


print(
    "\nNumber of connected components:",
    n_components,
)


print(
    "\nCurrent state order:"
)


for index, symbol in enumerate(
    analytic_sys.x_
):

    print(
        f"{index:2d}",
        f"{str(symbol):45s}",
        "component",
        labels[index],
    )


# --------------------------------------------------
# Display every mathematical component.
# --------------------------------------------------

print(
    "\n======================================"
)

print(
    "CONNECTED COMPONENTS"
)

print(
    "======================================"
)


for component in range(
    n_components
):

    indices = np.where(
        labels == component
    )[0]

    states = [
        str(
            analytic_sys.x_[index]
        )
        for index in indices
    ]

    contiguous = (
        len(indices) <= 1
        or np.all(
            np.diff(indices) == 1
        )
    )

    print(
        f"\nComponent {component}"
    )

    print(
        "indices:",
        list(indices),
    )

    print(
        "contiguous:",
        contiguous,
    )

    print(
        "states:"
    )

    for state in states:
        print(
            "   ",
            state,
        )


# --------------------------------------------------
# Test the CURRENT ODE-toolbox block routine.
# --------------------------------------------------

print(
    "\n======================================"
)

print(
    "CURRENT BLOCK-DIAGONAL ROUTINE"
)

print(
    "======================================"
)


try:

    blocks = (
        get_block_diagonal_blocks(
            A
        )
    )

    print(
        "\nSUCCESS"
    )

    print(
        "Number of blocks:",
        len(blocks),
    )

    print(
        "Block sizes:",
        [
            block.shape[0]
            for block in blocks
        ],
    )

except GetBlockDiagonalException:

    print(
        "\nFAILED"
    )

    print(
        "Current ordering causes "
        "GetBlockDiagonalException."
    )


# --------------------------------------------------
# Now REORDER states by connected component purely
# as a diagnostic experiment.
#
# This does NOT modify ODE-toolbox.
# --------------------------------------------------

permutation = np.argsort(
    labels,
    kind="stable",
)


A_reordered = A[
    np.ix_(
        permutation,
        permutation,
    )
]


print(
    "\n======================================"
)

print(
    "REORDERED COMPONENT TEST"
)

print(
    "======================================"
)


print(
    "\nReordered states:"
)


for new_index, old_index in enumerate(
    permutation
):

    print(
        f"{new_index:2d}",
        str(
            analytic_sys.x_[
                old_index
            ]
        ),
    )


try:

    reordered_blocks = (
        get_block_diagonal_blocks(
            A_reordered
        )
    )

    print(
        "\nREORDERED MATRIX: SUCCESS"
    )

    print(
        "Number of blocks:",
        len(reordered_blocks),
    )

    print(
        "Block sizes:",
        [
            block.shape[0]
            for block
            in reordered_blocks
        ],
    )

except GetBlockDiagonalException:

    print(
        "\nREORDERED MATRIX: "
        "STILL FAILED"
    )
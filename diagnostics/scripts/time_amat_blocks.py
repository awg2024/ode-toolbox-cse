import time

import numpy as np
import sympy

from amat_diagnostic_utils import (
    build_amat_system,
)

from odetoolbox.system_of_shapes import (
    get_block_diagonal_blocks,
)

from odetoolbox.config import Config


(
    model,
    shape_sys,
    analytic_sys,
    classification,
) = build_amat_system()


A = np.array(
    analytic_sys.A_,
    dtype=object,
)


# --------------------------------------------------
# Determine connected-component ordering
# --------------------------------------------------

connectivity = (
    (A != 0)
    |
    (A.T != 0)
)


import scipy.sparse


_, labels = (
    scipy.sparse.csgraph.connected_components(
        connectivity
    )
)


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


blocks = get_block_diagonal_blocks(
    A_reordered
)


h = sympy.Symbol(
    Config().output_timestep_symbol,
    real=True,
)


print(
    "Block sizes:",
    [
        block.shape[0]
        for block in blocks
    ],
)


# --------------------------------------------------
# Time each block separately
# --------------------------------------------------

total = 0.0


for index, block in enumerate(
    blocks
):

    block = sympy.Matrix(
        block
    )

    start = time.perf_counter()

    P_block = sympy.exp(
        block * h
    )

    elapsed = (
        time.perf_counter()
        - start
    )

    total += elapsed

    print(
        f"Block {index}",
        f"size={block.shape[0]}",
        f"time={elapsed:.3f}s",
    )


print(
    f"\nTotal block exponential time: "
    f"{total:.3f}s"
)
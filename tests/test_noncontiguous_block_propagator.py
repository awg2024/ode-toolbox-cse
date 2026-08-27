import sympy
from odetoolbox.system_of_shapes import (get_block_diagonal_blocks)
import numpy as np 
"""
Test propagator blcok equals propagator direct? 

"""



def test_noncontiguous_blocks_preserve_original_order():

    h = sympy.Symbol(
        "__h",
        real=True,
    )

    # Two independent 2x2 systems whose
    # states are interleaved:
    #
    # component A = indices 0, 2
    # component B = indices 1, 3
    #
    A = sympy.Matrix([ # use interleved matrix 
        [-1,  0,  1,  0],
        [ 0, -2,  0,  1],
        [ 0,  0, -1,  0],
        [ 0,  0,  0, -2],
    ])

    # Ground truth.
    expected = sympy.exp(A * h) # expected is the mathematical groundtruth

    blocks, permutation = (
        get_block_diagonal_blocks( # can you identify independent subsystems inside the matrix
            np.array(
                A,
                dtype=object,
            )
        )
    )

    block_propagators = [
        sympy.exp(
            sympy.Matrix(block)
            * h
        )
        for block in blocks
    ]

    reordered = sympy.diag(
        *block_propagators
    )

    inverse = np.argsort(
        permutation
    ).tolist()

    actual = reordered.extract(
        inverse,
        inverse,
    )

    for row in range(
        A.rows
    ):

        for col in range(
            A.cols
        ):

            assert sympy.simplify(
                expected[row, col]
                - actual[row, col]
            ) == 0
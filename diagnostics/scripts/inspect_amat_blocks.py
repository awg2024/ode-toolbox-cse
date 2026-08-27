import numpy as np
import scipy.sparse

from amat_diagnostic_utils import (
    build_amat_system,
)

from odetoolbox.system_of_shapes import (
    get_block_diagonal_blocks,
    GetBlockDiagonalException,
)


"""
Here we select a specific analytical subsystem A and we specifically
investigating why block optimisation fails in odetoolbox. 

It builds a complete diagnostic simulation that maps out the exact variables of your AMAT model, 
tests if the current odetoolbox code crashes on it, and then conducts a reordering experiment to show
how a simple matrix permutation makes the optimization work.
"""

# calling helper function 
(model,shape_sys,analytic_sys,classification,) = build_amat_system()

print("Only selecting matrix A:", analytic_sys.A_)
# Work specifically with the analytical subsystem A.
A = np.array(analytic_sys.A_,dtype=object)
print("analytical system matrix") # shape of the system, given there are 8 state variables in amat ---> 8x8 matrix 
print("Matrix A shape:", A.shape)

# Construct the exact undirected connectivity graph
# used by get_block_diagonal_blocks().
connectivity = ((A != 0) | (A.T != 0)) # checks which variables talk to eachother binary 
n_components, labels = (scipy.sparse.csgraph.connected_components(connectivity)) # split them into cluster islands 

# ignore depedencies (direction) 
print("\nNumber of connected components:",n_components)
print("\nCurrent state order:")

for index, symbol in enumerate(analytic_sys.x_):

    print(f"{index:2d}",
        f"{str(symbol):45s}",
        "component",
        labels[index])

# Display every mathematical component of this amat matrix A
# Here we are displaying the full 8x8 matrix without any singularity branching occuring yet 
print("connected components")

for component in range(n_components):

    
    indices = np.where(labels == component)[0]

    states = [str(analytic_sys.x_[index]) for index in indices]

    # calculate the stat variables belonging to those islands 
    # if they are together they are contiguous: true or else contiguous: false 
    contiguous = (len(indices) <= 1 or np.all(np.diff(indices) == 1))

    print(f"\nComponent {component}")
    print("indices:", list(indices))

    print("are they contiguous?",contiguous,)

    print("states:")
    for state in states:
        print("   ",state)



print("current block-diagonal route, with checks for contiguous matrices")

try:

    blocks = (get_block_diagonal_blocks(A))

    print("success in block diagonal matrix is successfully contiguous")
    print("Number of blocks:",len(blocks))
    print("Block sizes:",[block.shape[0] for block in blocks])

except GetBlockDiagonalException:

    print("fail in block diagonal matrix is not contiguous")
    print(
        "Current ordering causes "
        "GetBlockDiagonalException."
    )

# modifying and reordering states by connected components to make 
# it contiguous (diagnostic experiment to see if we can get this function to pass)


# sort stable variables based on labels 
permutation = np.argsort(labels,kind="stable")

# reorder the state variable
A_reordered = A[np.ix_(permutation,permutation)]


print(
    "REORDERED COMPONENT TEST"
)

print(
    "======================================"
)

print("attempt 2: current block-diagonal route, with checks for contiguous matrices")
print("reordered states")

for new_index, old_index in enumerate(permutation):

    print(f"{new_index:2d}",str(analytic_sys.x_[old_index]))

try:

    reordered_blocks = (
        get_block_diagonal_blocks(
            A_reordered
        )
    )

    print("success in block diagonal matrix is successfully contiguous")
    print("Number of reordered blocks:",len(reordered_blocks))
    print("Block sizes:",[block.shape[0] for block in reordered_blocks])
    
except GetBlockDiagonalException:

    print("fail in block diagonal matrix is not contiguous")
    print(
        "Current ordering causes "
        "GetBlockDiagonalException."
    )
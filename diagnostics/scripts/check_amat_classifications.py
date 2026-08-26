from amat_diagnostic_utils import (
    build_amat_system,
)


(
    model,
    shape_sys,
    analytic_sys,
    classification,
) = build_amat_system()


print(
    "\n======================================"
)

print(
    "AMAT ANALYTICAL CLASSIFICATION"
)

print(
    "======================================\n"
)


for symbol, solvable in (
    classification.items()
):

    status = (
        "ANALYTICAL"
        if solvable
        else "NUMERICAL"
    )

    print(
        f"{str(symbol):45s} {status}"
    )


print(
    "\nAnalytical subsystem size:",
    len(analytic_sys.x_),
)


print(
    "\nAnalytical states:"
)

for symbol in analytic_sys.x_:
    print(
        "   ",
        symbol,
    )
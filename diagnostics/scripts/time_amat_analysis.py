import argparse
import json
import time

import odetoolbox


parser = argparse.ArgumentParser()

parser.add_argument(
    "mode",
    choices=[
        "base",
        "detect",
        "full",
        "alt-base",
    ],
)

args = parser.parse_args()


with open(
    "amat_json/amat_ode_toolbox_input.json",
    "r",
    encoding="utf-8",
) as file:

    model = json.load(file)


kwargs = {
    "disable_stiffness_check": True,

    # Absolutely no CSE during diagnosis.
    "enable_cse": False,
}


if args.mode == "base":

    print(
        "\nMODE: BASE PROPAGATOR ONLY"
    )

    kwargs[
        "disable_singularity_detection"
    ] = True


elif args.mode == "detect":

    print(
        "\nMODE: PROPAGATOR + "
        "SINGULARITY DETECTION"
    )

    kwargs[
        "disable_singularity_detection"
    ] = False

    kwargs[
        "disable_singularity_mitigation"
    ] = True


elif args.mode == "full":

    print(
        "\nMODE: FULL SINGULARITY "
        "MITIGATION"
    )

    kwargs[
        "disable_singularity_detection"
    ] = False

    kwargs[
        "disable_singularity_mitigation"
    ] = False


elif args.mode == "alt-base":

    print(
        "\nMODE: ALTERNATIVE MATRIX "
        "EXPONENTIAL"
    )

    kwargs[
        "disable_singularity_detection"
    ] = True

    kwargs[
        "use_alternative_expM"
    ] = True


start = time.perf_counter()


result = odetoolbox.analysis(
    model,
    **kwargs,
)


elapsed = (
    time.perf_counter()
    - start
)


print(
    f"\nElapsed: {elapsed:.3f} seconds"
)


print(
    "\nReturned solver blocks:",
    len(result),
)


for index, solver in enumerate(
    result
):

    print(
        f"\nSolver {index}:",
        solver["solver"],
    )

    print(
        "States:",
        solver[
            "state_variables"
        ],
    )

    if "conditions" in solver:

        print(
            "Final merged conditions:",
            len(
                solver[
                    "conditions"
                ]
            ),
        )
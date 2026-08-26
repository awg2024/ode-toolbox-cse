import json
import sympy

from odetoolbox import (
    _read_global_config,
    _from_json_to_shapes,
    _find_analytically_solvable_equations,
)

from odetoolbox.system_of_shapes import (
    SystemOfShapes,
)


INPUT_FILE = (
    "scripts/tests/"
    "amat_ode_toolbox_input.json"
)


with open(
    INPUT_FILE,
    "r",
    encoding="utf-8",
) as file:

    model = json.load(file)


# -----------------------------------------
# Process global configuration
# -----------------------------------------

_read_global_config(
    model
)


# -----------------------------------------
# Convert parameters to SymPy symbols
# -----------------------------------------

parameters = {}

for name, value in model[
    "parameters"
].items():

    parameters[
        sympy.Symbol(
            name,
            real=True,
        )
    ] = value


# -----------------------------------------
# Build shapes and complete state system
# -----------------------------------------

shapes, parameters = (
    _from_json_to_shapes(
        model,
        parameters=parameters,
    )
)


shape_sys = (
    SystemOfShapes.from_shapes(
        shapes,
        parameters=parameters,
    )
)


# -----------------------------------------
# ONLY perform analytical classification
# -----------------------------------------

dependency_edges, classification = (
    _find_analytically_solvable_equations(
        shape_sys,
        shapes,
        parameters=parameters,
    )
)


print("\n================================")
print("DEPENDENCY EDGES")
print("================================\n")

print(dependency_edges)

print("\n================================")
print("CLASSIFICATION")
print("================================\n")

for symbol, solvable in classification.items():
    print(
        f"{str(symbol):35s}",
        "ANALYTICAL" if solvable else "NUMERICAL",
    )

print("\n================================")
print("DEPENDENCY GRAPH")
print("================================\n")

for src, dst in dependency_edges:
    print(f"{src}  --->  {dst}")

print("\n================================")
print("INCOMING DEPENDENCIES")
print("================================\n")

nodes = set()
for src, dst in dependency_edges:
    nodes.add(src)
    nodes.add(dst)

for node in nodes:
    incoming = [
        src
        for src, dst in dependency_edges
        if dst == node and src != dst
    ]

    outgoing = [
        dst
        for src, dst in dependency_edges
        if src == node and src != dst
    ]

    print(f"{node}")
    print(f"    incoming: {incoming}")
    print(f"    outgoing: {outgoing}")



print("\n================================")
print("DEPENDENCY SUMMARY")
print("================================\n")

for node in sorted(nodes, key=str):
    incoming = sorted(
        {
            src
            for src, dst in dependency_edges
            if dst == node and src != dst
        },
        key=str,
    )

    outgoing = sorted(
        {
            dst
            for src, dst in dependency_edges
            if src == node and src != dst
        },
        key=str,
    )

    print(f"{node}")
    print(f"    # incoming: {len(incoming)}")
    print(f"    incoming:  {incoming}")
    print(f"    # outgoing: {len(outgoing)}")
    print(f"    outgoing:  {outgoing}")
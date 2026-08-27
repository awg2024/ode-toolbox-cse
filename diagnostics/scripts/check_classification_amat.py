import json
import sympy
from odetoolbox import (
    _read_global_config,
    _from_json_to_shapes,
    _find_analytically_solvable_equations)

from odetoolbox.system_of_shapes import (SystemOfShapes)

# input json file, 
INPUT_FILE = ("json_testfile/amat_ode_toolbox_input.json")


with open(INPUT_FILE, "r", encoding="utf-8") as file:
    model = json.load(file)

# Process global configuration from __init__ looking for "options" in keys 
_read_global_config(model)

# Convert parameters to SymPy symbols matching spec and enabling model run
parameters = {}

for name, value in model["parameters"].items():
    parameters[sympy.Symbol(name,real=True)] = value

# Converting json to shapes / state variables?  
shapes, parameters = (
    _from_json_to_shapes(
        model,
        parameters=parameters,
    )
)

# Construct the global system matrix
shape_sys = (SystemOfShapes.from_shapes(shapes, parameters=parameters))

# ONLY perform analytical classification (we know this is the case since it's an amat model)
dependency_edges, classification = (
    _find_analytically_solvable_equations(
        shape_sys,
        shapes,
        parameters=parameters,
    )
)

print("Dependency edges or the mathematical causal links between different state variables:")
print(dependency_edges)

print("classification of how these state variables were solved:")
for symbol, solvable in classification.items():
    print(
        f"{str(symbol):35s}",
        "ANALYTICAL" if solvable else "NUMERICAL",
    )

print("---")
print("src (Source): The independent variable or the shape that provides the data.")
print("dst (Destination): The dependent variable or the shape that relies on src.")
for src, dst in dependency_edges:
    print(f"{src}  --->  {dst}")
print("---")

print("incoming src dependencies: ")
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

    print("---")
    print(f"{node}")
    print(f"    incoming: {incoming}")
    print(f"    outgoing: {outgoing}")

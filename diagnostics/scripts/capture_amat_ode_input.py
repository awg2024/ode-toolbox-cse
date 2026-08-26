from pathlib import Path
from pynestml.codegeneration.nest_code_generator_utils import NESTCodeGeneratorUtils


NESTML_MODEL = Path("../tests/amat_neuron.nestml") # hard code path 

if not NESTML_MODEL.exists():
    raise FileNotFoundError(NESTML_MODEL.resolve())

print("processing:", NESTML_MODEL.resolve())

NESTCodeGeneratorUtils.generate_code_for(str(NESTML_MODEL), module_name="amat_capture_module", logging_level="INFO") # load amat in cse odetoolbox env 
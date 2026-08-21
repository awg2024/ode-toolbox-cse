#
# test_disabled_preserves_legacy_cse.py
#
# This file is part of the NEST ODE toolbox.
#
# Copyright (C) 2017 The NEST Initiative
#
# The NEST ODE toolbox is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 2 of
# the License, or (at your option) any later version.
#
# The NEST ODE toolbox is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Existing users upstream should notice no difference in their output. 
"""


import pytest
import odetoolbox
from tests.test_utils import load_test_json

def test_cse_disabled_preserves_legacy():
    
     # mocking a .json from a model description file 
    indict = {
        "dynamics": [ 
            {
                "expression":
                    "x' = -x / tau", 
                    "initial_value":
                        "1",
            }   
        ]    
    }

    default_result = odetoolbox.analysis( # default cse is off by default 
        indict,
        disable_stiffness_check=True,
        disable_singularity_detection=True, 
    )

    explicitly_disabled_result = odetoolbox.analysis( # explicitly false cse 
        indict,
        disable_stiffness_check=True,
        disable_singularity_detection=True,
        enable_cse=False
    )

    assert(default_result == explicitly_disabled_result)
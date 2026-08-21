#
# test_numeric_solver_cse.py
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

""""
This script provides a tests the solving ability of linear numerical aeif_cond_beta ODEs with CSE optimisation 
"""

import sympy 
from odetoolbox.expression_optimisation import apply_cse_to_solver


def test_numeric_solver_cse():

    solver_dict = {
        {
            "solver": "analytical",  # Block 0: Tool tells us it is analytical
            "state_variables": ["g_in__X__inh_spikes", "g_in__DOLLAR__X__inh_spikes", "g_ex__X__exc_spikes", "g_ex__DOLLAR__X__exc_spikes"],
            "propagators": {
                "__P__g_ex__DOLLAR__X__exc_spikes__g_ex__DOLLAR__X__exc_spikes": "exp(-__h/tau_syn_decay_E)",
                "__P__g_ex__X__exc_spikes__g_ex__DOLLAR__X__exc_spikes": "-tau_syn_decay_E*tau_syn_rise_E*exp(-__h/tau_syn_rise_E)/(tau_syn_decay_E - tau_syn_rise_E) + tau_syn_decay_E*tau_syn_rise_E*exp(-__h/tau_syn_decay_E)/(tau_syn_decay_E - tau_syn_rise_E)",
                "__P__g_ex__X__exc_spikes__g_ex__X__exc_spikes": "exp(-__h/tau_syn_rise_E)",
                "__P__g_in__DOLLAR__X__inh_spikes__g_in__DOLLAR__X__inh_spikes": "exp(-__h/tau_syn_decay_I)",
                "__P__g_in__X__inh_spikes__g_in__DOLLAR__X__inh_spikes": "-tau_syn_decay_I*tau_syn_rise_I*exp(-__h/tau_syn_rise_I)/(tau_syn_decay_I - tau_syn_rise_I) + tau_syn_decay_I*tau_syn_rise_I*exp(-__h/tau_syn_decay_I)/(tau_syn_decay_I - tau_syn_rise_I)",
                "__P__g_in__X__inh_spikes__g_in__X__inh_spikes": "exp(-__h/tau_syn_rise_I)"
            },
            "update_expressions": {
                "g_ex__DOLLAR__X__exc_spikes": "__P__g_ex__DOLLAR__X__exc_spikes__g_ex__DOLLAR__X__exc_spikes*g_ex__DOLLAR__X__exc_spikes",
                "g_ex__X__exc_spikes": "__P__g_ex__X__exc_spikes__g_ex__DOLLAR__X__exc_spikes*g_ex__DOLLAR__X__exc_spikes + __P__g_ex__X__exc_spikes__g_ex__X__exc_spikes*g_ex__X__exc_spikes",
                "g_in__DOLLAR__X__inh_spikes": "__P__g_in__DOLLAR__X__inh_spikes__g_in__DOLLAR__X__inh_spikes*g_in__DOLLAR__X__inh_spikes",
                "g_in__X__inh_spikes": "__P__g_in__X__inh_spikes__g_in__DOLLAR__X__inh_spikes*g_in__DOLLAR__X__inh_spikes + __P__g_in__X__inh_spikes__g_in__X__inh_spikes*g_in__X__inh_spikes"
            }
        },
        {
            "solver": "numeric",  # <-- Block 1: Tool tells us it is numeric
            "state_variables": ["V_m", "w"],
            "update_expressions": {
                "V_m": "((-g_L) * ((min(V_m, V_peak)) - E_L) + (g_L * Delta_T * exp((((min(V_m, V_peak)) - V_th) / Delta_T))) - (g_ex__X__exc_spikes * 1.0 * ((min(V_m, V_peak)) - E_exc)) - (g_in__X__inh_spikes * 1.0 * ((min(V_m, V_peak)) - E_inh)) - w + I_e + I_stim) / C_m",
                "w": "(a * ((min(V_m, V_peak)) - E_L) - w) / tau_w"
            }
        }
    }
    

    results = []
    result = apply_cse_to_solver(solver_dict)

    print(f"\nAfter CSE Output Dictionary: {result}") 
  
    assert result["solver"] == "numeric", "Failed: The solver type mutated or was lost."
    assert "cse" in result, "Failed: 'cse' key dictionary was never initialized."

    assert "update_expressions" in result["cse"], (
        "Failed: 'update_expressions' was not processed or missing inside the inner cse tracker.")

    assert "update_expressions" in result, "Failed: Top-level update expressions dictionary was stripped."

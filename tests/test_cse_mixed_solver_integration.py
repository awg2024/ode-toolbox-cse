# test_cse_mixed_solver_integration.py
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

from tests.test_utils import load_test_json

from tests.cse_test_utils import (
    assert_cse_dependency_order,
    assert_cse_region_equivalent,
    assert_cse_region_profitable,
    assert_cse_serialized,
    assert_cse_temporaries_disjoint,
    assert_shape_structure_preserved,
    assert_solver_metadata_preserved,
    get_solver,
    run_cse_analysis_pair,
)


class TestCSEMixedSolver:
    """
    Isolated ODE-toolbox validation of CSE for a system containing both an
    analytical solver block and a numerical solver block.
    """

    def test_mixed_cse_pipeline(self):

        indict = load_test_json("cse_mixed.json")

        pair = run_cse_analysis_pair(
            indict,
            disable_stiffness_check=True,
            disable_singularity_detection=True,
            log_level="DEBUG",
        )

    
        # Mixed system must produce exactly two solver blocks.
        assert len(pair.baseline_solvers) == 2
        assert len(pair.cse_solvers) == 2

    
        # Extract corresponding blocks by solver type rather than relying on
        # list ordering.
        baseline_analytical = get_solver(
            pair.baseline_solvers,
            "analytical")

        baseline_numerical = get_solver(
            pair.baseline_solvers,
            "numeric")

        cse_analytical = get_solver(
            pair.cse_solvers,
            "analytical")

        cse_numerical = get_solver(
            pair.cse_solvers,
            "numeric")

        # Baseline should have no CSE metadata.
        assert ("cse" not in baseline_analytical)
        assert ("cse" not in baseline_numerical)

        # Internal ODE representation cannot change.
        assert_shape_structure_preserved(pair)

        # Solver partition / metadata cannot chang
        assert_solver_metadata_preserved(
            baseline_analytical,
            cse_analytical,
        )

        assert_solver_metadata_preserved(
            baseline_numerical,
            cse_numerical,
        )

    
        # copied in from analytical test checks ;;;
        
        assert "cse" in cse_analytical

        assert (
            "propagators"
            in cse_analytical["cse"]
        )

        assert_cse_region_equivalent(
            baseline_analytical,
            cse_analytical,
            "propagators",
            require_cse=True,
        )

        assert_cse_region_profitable(
            baseline_analytical,
            cse_analytical,
            "propagators",
        )

        # copied in from numerical test checks ;;; 

        assert "cse" in cse_numerical

        assert (
            "update_expressions"
            in cse_numerical["cse"]
        )

        assert_cse_region_equivalent(
            baseline_numerical,
            cse_numerical,
            "update_expressions",
            require_cse=True,
        )

        assert_cse_region_profitable(
            baseline_numerical,
            cse_numerical,
            "update_expressions",
        )

        
        # Both solver blocks must have had their metadata serialized.
        # THIS directly catchesthe current serialization-loop problem.
        assert_cse_serialized(cse_analytical)
        assert_cse_serialized(cse_numerical)

        
        # CSE namespace isolation, ensure tmp variables are not leaking or colliding with other branches (e.g., numerical vs analytical)
        assert_cse_temporaries_disjoint(pair.cse_solvers)

        # Each solver's CSE dependency order is valid.
        assert_cse_dependency_order(cse_analytical)
        assert_cse_dependency_order(cse_numerical)
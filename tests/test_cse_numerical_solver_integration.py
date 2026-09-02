#
# test_cse_numerical_solver_integration.py
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
    assert_shape_structure_preserved,
    assert_solver_metadata_preserved,
    get_solver,
    run_cse_analysis_pair,
)


class TestCSENumericalSolver:
    """
    Isolated ODE-toolbox validation of CSE applied to a numerical solver.
    """

    def test_numerical_cse_pipeline(self):

        indict = load_test_json("cse_numerical.json")

        pair = run_cse_analysis_pair(
            indict,
            disable_stiffness_check=True,
            disable_singularity_detection=True,
            log_level="DEBUG",
        )

        # Nonlinear json should produce only a numerical solver.
        assert len(pair.baseline_solvers) == 1
        assert len(pair.cse_solvers) == 1

        baseline_solver = get_solver(pair.baseline_solvers,"numeric")
        cse_solver = get_solver(pair.cse_solvers,"numeric")

        assert "cse" not in baseline_solver

    
        # CSE must not affect ODE parsing / Shape construction.
        assert_shape_structure_preserved(pair)

        # Solver selection, state variables, etc. must remain unchanged.
        assert_solver_metadata_preserved(baseline_solver, cse_solver)


        # Numerical solver has no propagators.
        assert "propagators" not in baseline_solver
        assert "propagators" not in cse_solver

        # Numerical CSE operates on RHS/update expressions.
        assert "update_expressions" in baseline_solver
        assert "update_expressions" in cse_solver

        # cse is present since we delibately made repeated expressions
        assert "cse" in cse_solver
        assert ("update_expressions" in cse_solver["cse"])

        # Reconstruct CSE expressions and prove equivalence against the
        # unoptimised numerical RHS.
        assert_cse_region_equivalent(
            baseline_solver,
            cse_solver,
            "update_expressions",
            require_cse=True)

        # Transformation must be genuinely cheaper.
        assert_cse_region_profitable(baseline_solver, cse_solver, "update_expressions")

        #  Check that the metadata returned by _analysis must is JSON-safe.
        assert_cse_serialized(cse_solver)

        # Temporary dependency order must be valid (e.g., not calling a tmpvar before it's defined)
        assert_cse_dependency_order(cse_solver)

    def test_legacy_iaf_cond_alpha_preserved(self):

        indict = load_test_json("iaf_cond_alpha.json")

        pair = run_cse_analysis_pair(
            indict,
            disable_stiffness_check=True,
            disable_analytic_solver=True,
            disable_singularity_detection=True,
            log_level="DEBUG")
            
        baseline_solver = get_solver(
            pair.baseline_solvers,
            "numeric")

        cse_solver = get_solver(
            pair.cse_solvers,
            "numeric")

        assert_solver_metadata_preserved(
            baseline_solver,
            cse_solver)

        assert_cse_region_equivalent(
            baseline_solver, # check these two solvers are equivalant and upstream numerical users are not being affected. 
            cse_solver,
            "update_expressions",
            require_cse=False)


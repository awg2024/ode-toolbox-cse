#
# test_cse_analytic_solver_integration.py
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
    run_cse_analysis_pair)


class TestCSEAnalyticalSolver:
    """
    Isolated ODE-toolbox validation of analytical CSE.

    The same model is analysed with CSE disabled and enabled. The test verifies
    that CSE changes the symbolic representation but not its mathematical
    meaning.
    """

    def test_analytical_cse_pipeline(self):

        indict = load_test_json("cse_analytical.json")

        pair = run_cse_analysis_pair(
            indict,
            disable_stiffness_check=True,
            disable_singularity_detection=True,
            log_level="DEBUG",
        )

        # ensure there is output for both conditions 
        assert len(pair.baseline_solvers) == 1
        assert len(pair.cse_solvers) == 1

        # ensure solver specified is analytical 
        baseline_solver = get_solver(pair.baseline_solvers, "analytical")
        cse_solver = get_solver(pair.cse_solvers,"analytical")

        # CSE should not exist in baseline.
        assert "cse" not in baseline_solver


        # Internal Shape/SystemOfShapes construction must be identical.
        assert_shape_structure_preserved(pair)

        # Solver metadata cannot change.
        assert_solver_metadata_preserved(baseline_solver, cse_solver)

        # json file has repeated propagators (==> profitable) so ensure cse has occured
        assert "cse" in cse_solver

        assert ("propagators" in cse_solver["cse"])

        # reconstruct every CSE expression and prove that it is mathematically the same to the non-cse propagator
        assert_cse_region_equivalent(baseline_solver, cse_solver, "propagators", require_cse=True)

        # verifies that the optimiser acted correctly, reducing mathematical operation symbol counts. 
        assert_cse_region_profitable(baseline_solver, cse_solver, "propagators")

        # Check that the metadata returned by _analysis must is JSON-safe.
        assert_cse_serialized(cse_solver)

        # ensure that CSE solver doesnt contain forward references (e.g., calling a tmp_var before it's defined)
        assert_cse_dependency_order(cse_solver)
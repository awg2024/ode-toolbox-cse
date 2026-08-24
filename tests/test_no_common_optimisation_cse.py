"""
This test investigates what happens when you apply cse to an already optimised solver. 
"""




# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.
#

""""
This script provides a test for an analytical solver cse, utilising the generate_propagator_solver()



"""

import sympy 
from odetoolbox.expression_optimisation import apply_cse_to_solver, restore_cse_expression
from .context import odetoolbox

def test_cse_no_common_expression():
    x, y = sympy.symbols("x y", real=True)

    indict = {
        "dynamics": [ # mocking indict that will be passed into _analysis 
            { 
                "expression": "x' = x + 1",
                "initial_value": "1",
            },
            {
                "expression": "y' = y + 1", 
                "initial_value": "2",
            }
        ]
    }


    # Explicitly testing the False behavior
    result = odetoolbox.analysis(
        indict, 
        enable_cse=True
    )


    # Check that 'cse' key was not added to the root dictionary
    assert "cse" not in result 
    
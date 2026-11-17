import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "src"))
import numpy as np

from parse_case import parse_case, bus_index_map
from ybus import build_ybus
from run_veragrid import veragrid_ybus

def test_ybus():
    case = parse_case("data/case9.m")

    my_ybus = build_ybus(case)
    vg_ybus, vg_bus_numbers = veragrid_ybus("data/case9.m")

    my_bus_index = bus_index_map(case['bus'])

    order = [
        my_bus_index[bus] for bus in vg_bus_numbers
    ]

    ybus_organized = my_ybus[np.ix_(order, order)]

    max_diff = np.max(np.abs(ybus_organized - vg_ybus))

    tol = 1e-8 # absolute error for p.u. admittance; rounding is about 1e-16 for 64-bit floats; smallest real modelling error (dropping jb/2) is about 0.08"

    print(max_diff)

    assert max_diff < tol, (
        f"Ybus matrices differ by {max_diff}, which is greater than the tolerance of {tol}"
    )
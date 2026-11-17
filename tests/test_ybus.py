"""Check my Ybus against VeraGrid's on case9.

What a pass shows: the branch pi models and bus alignment match VeraGrid.
What it does not show: bus shunts, out-of-service branches and taps, because
case9 has none of them.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "src")) # let Python find src/
import numpy as np

from parse_case import parse_case, bus_index_map
from ybus import build_ybus
from run_veragrid import veragrid_ybus

def test_ybus():
    case = parse_case("data/case9.m")

    my_ybus = build_ybus(case)
    vg_ybus, vg_bus_numbers = veragrid_ybus("data/case9.m")

    my_bus_index = bus_index_map(case['bus']) # bus number -> my row

    order = [
        my_bus_index[bus] for bus in vg_bus_numbers # my row for each of VeraGrid's buses, in VeraGrid's order
    ]

    ybus_organized = my_ybus[np.ix_(order, order)] # reorder my rows and columns to match VeraGrid

    max_diff = np.max(np.abs(ybus_organized - vg_ybus)) # largest difference over all entries, p.u.

    # absolute error in p.u. admittance. Floating-point rounding is ~1e-16 relative,
    # ~1e-14 for entries of 10-40 (measured: 4.4e-16). The smallest real modelling
    # error (dropping j*b/2 on one line) is ~0.08. 1e-8 sits well between the two.
    tol = 1e-8

    print(max_diff) # shown with: python -m pytest -v -s

    assert max_diff < tol, (
        f"Ybus matrices differ by {max_diff}, which is greater than the tolerance of {tol}"
    )

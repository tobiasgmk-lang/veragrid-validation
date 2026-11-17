import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "src"))
import numpy as np
import pytest

from parse_case import parse_case, BUS_I
from ybus import build_ybus
from run_veragrid import veragrid_ybus

def test_ybus():
    case = parse_case("data/case9.m")

    my_ybus = build_ybus(case)
    vg_ybus, vg_bus_numbers = veragrid_ybus("data/case9.m")

    my_bus_numbers = case['bus'][:, BUS_I].astype(int)

    my_bus_index = {}

    for i in range(len(my_bus_numbers)):
        bus_number = int(my_bus_numbers[i])
        my_bus_index[bus_number] = i

    order = [
        my_bus_index[bus] for bus in vg_bus_numbers
    ]

    ybus_organized = my_ybus[np.ix_(order, order)]

    max_diff = np.max(np.abs(ybus_organized - vg_ybus))

    tol = 1e-8 # tolerance, 6 orders of mag above nose, 7 below a real error

    print(max_diff)

    assert max_diff < tol, (
        f"Ybus matrices differ by {max_diff}, which is greater than the tolerance of {tol}"
    )
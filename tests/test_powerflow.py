"""Compare my Newton-Raphson power flow with VeraGrid's on case9.

Both solvers read the same MATPOWER file, start from the same flat start and stop
when the largest power mismatch is below 1e-6 p.u. Any difference beyond that
noise comes from the method or a convention, not the data.
Results are written to results/case9_comparison.csv and results/case9_summary.csv.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "src")) # let Python find src/
import numpy as np

from parse_case import parse_case, bus_index_map
from newton_raphson import solve_power_flow
from run_veragrid import run_veragrid
from report import write_comparison

def test_powerflow():
    case = parse_case("data/case9.m")

    # run both solvers on the same case
    mine = solve_power_flow(case)
    vg = run_veragrid("data/case9.m")

    # both solvers should converge
    assert mine["converged"], "mine did not converge"
    assert vg["converged"], "veragrid did not converge"

    # align the bus results by bus number
    my_bus_index = bus_index_map(case['bus'])

    order = [
            my_bus_index[bus] for bus in vg["bus"] # my row for each of VeraGrid's buses, in VeraGrid's order
        ]

    # reorder results to match VeraGrid's bus order
    final_Vm = mine["Vm"][order]
    final_Va = mine["Va_deg"][order]

    # compare voltage magnitudes and angles
    max_dVm = np.max(np.abs(final_Vm - vg["Vm"]))
    delta_A = (final_Va - vg["Va_deg"] + 180) % 360 - 180 # wrap to [-180, 180]
    max_dVa = np.max(np.abs(delta_A))

    # compare active power losses
    delta_losses = abs(mine["losses_MW"] - vg["losses_MW"])

    # tolerances for comparison
    # each tolerance sits between two numbers:
    #   noise floor: both solvers stop once the largest power mismatch < 1e-6 p.u., so each
    #                answer can sit ~1e-6 away from the true solution even if both are correct
    #   real error:  the smallest modelling mistake worth catching, e.g. dropping j*b/2 on one
    #                line, moves |V| by ~1e-3 p.u. (estimated, not measured)
    vm_tol = 1e-5 # p.u.: noise ~1e-6 p.u., so 10x above noise and ~100x below a real error
    va_tol = 1e-3 # degrees: noise ~1e-6 rad = ~6e-5 deg, so ~17x above noise
    loss_tol = 1e-2 # MW: noise = 8 buses x 1e-6 p.u. x 100 MVA = 8e-4 MW, so ~12x above noise

    print(f"Maximum |∆V|: {max_dVm:.2e} p.u.")
    print(f"Maximum |∆θ|: {max_dVa:.2e} degrees")
    print(f"|∆losses|: {delta_losses:.2e} MW")
    print(f"My iterations: {mine['iterations']}, VeraGrid iterations: {vg['iterations']}")

    # write the CSVs before asserting, so a failing run still leaves results to inspect
    write_comparison("case9", vg["bus"], final_Vm, vg["Vm"], final_Va, vg["Va_deg"], mine, vg)

    # assert voltage within tolerances
    assert max_dVm < vm_tol, (
        f"Voltage magnitude difference {max_dVm:.2e} exceeds tolerance {vm_tol:.2e}"
    )
    assert max_dVa < va_tol, (
        f"Voltage angle difference {max_dVa:.2e} exceeds tolerance {va_tol:.2e}"
    )
    assert delta_losses < loss_tol, (
        f"Losses difference {delta_losses:.2e} exceeds tolerance {loss_tol:.2e}"
    )
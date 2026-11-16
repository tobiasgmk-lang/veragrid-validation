"""Run VeraGrid's AC power flow on a MATPOWER case and report results by bus number."""
import sys

import numpy as np
import VeraGridEngine.api as gce


def run_veragrid(case_path):
    grid = gce.open_file(case_path)
    options = gce.PowerFlowOptions(
        solver_type=gce.SolverType.NR,   # Newton-Raphson, same method as mine
        tolerance=1e-6,                  # VeraGrid default, stated explicitly
        max_iter=25,
        control_q=False,                 # reactive limits OFF for a like-for-like comparison
        retry_with_other_methods=False,  # if NR fails, fail loudly instead of switching solver
    )
    res = gce.power_flow(grid, options)

    V = np.asarray(res.voltage)                     # complex bus voltages, p.u.
    bus_numbers = [int(b.name) for b in grid.buses]  # MATPOWER bus number for each row
    return {
        "bus": bus_numbers,
        "Vm": np.abs(V),
        "Va_deg": np.degrees(np.angle(V)),
        "converged": bool(res.converged),
        "iterations": res.iterations,
        "error": res.error,
        "losses_MW": float(np.sum(np.real(res.losses))),
    }


def veragrid_ybus(case_path):
    """VeraGrid's own Ybus (dense complex array) and the bus number of each row."""
    grid = gce.open_file(case_path)
    nc = gce.compile_numerical_circuit_at(grid)      # VeraGrid's internal numeric model
    Ybus = nc.get_admittance_matrices().Ybus.toarray()  # sparse -> ordinary NumPy array
    bus_numbers = [int(b.name) for b in grid.buses]
    return Ybus, bus_numbers


if __name__ == "__main__":
    out = run_veragrid(sys.argv[1])
    print(f"converged={out['converged']}  iterations={out['iterations']}  "
          f"error={out['error']:.2e}  losses={out['losses_MW']:.3f} MW")
    print(" bus   |V| (p.u.)   angle (deg)")
    for b, vm, va in sorted(zip(out["bus"], out["Vm"], out["Va_deg"])):
        print(f"{b:4d}   {vm:9.4f}   {va:10.4f}")

"""Write comparison results to CSV files in results/.

Two files per case:
    results/<case>_comparison.csv   one row per bus: both solvers' |V| and angle, and the differences
    results/<case>_summary.csv      one row per solver: converged, iterations, final mismatch, losses
"""
import csv
import pathlib

import numpy as np

RESULTS = pathlib.Path(__file__).parents[1] / "results" # results/ next to src/

def write_comparison(case_name, bus, Vm_mine, Vm_vg, Va_mine, Va_vg, mine, vg):
    """bus, Vm_*, Va_* must already be aligned by bus number (same order for both solvers).

    mine and vg are the dicts returned by solve_power_flow and run_veragrid.
    """
    RESULTS.mkdir(exist_ok=True)

    with open(RESULTS / f"{case_name}_comparison.csv", "w", newline="") as f: # "w" overwrites the file every run
        w = csv.writer(f)
        # 8 decimals for values so tiny differences stay visible; differences in scientific notation
        w.writerow(["bus", "Vm_mine_pu", "Vm_veragrid_pu", "dVm_pu",
                    "Va_mine_deg", "Va_veragrid_deg", "dVa_deg"])
        for b, vm1, vm2, va1, va2 in zip(bus, Vm_mine, Vm_vg, Va_mine, Va_vg):
            w.writerow([b, f"{vm1:.8f}", f"{vm2:.8f}", f"{abs(vm1 - vm2):.3e}",
                        f"{va1:.8f}", f"{va2:.8f}", f"{abs(va1 - va2):.3e}"])

    with open(RESULTS / f"{case_name}_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["solver", "converged", "iterations", "final_max_mismatch_pu",
                    "losses_MW", "max_dVm_pu", "max_dVa_deg"])
        max_dV = np.max(np.abs(np.asarray(Vm_mine) - np.asarray(Vm_vg)))
        max_dA = np.max(np.abs(np.asarray(Va_mine) - np.asarray(Va_vg)))
        for name, r in (("mine", mine), ("veragrid", vg)):
            w.writerow([name, r["converged"], r["iterations"], f"{r['error']:.3e}",
                        f"{r['losses_MW']:.6f}", f"{max_dV:.3e}", f"{max_dA:.3e}"])

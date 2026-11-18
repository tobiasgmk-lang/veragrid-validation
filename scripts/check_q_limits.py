"""Step 4 check: do generator reactive (Q) limits change VeraGrid's case9 answer?

Runs VeraGrid twice (Q limits off, then on), prints the largest change in |V| and
angle, and prints each generator's Q against its limits so the margin is visible.
Run from the project root:  python scripts/check_q_limits.py
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "src")) # let Python find src/
import numpy as np

from parse_case import parse_case, bus_index_map, GEN_BUS, QMAX, QMIN, QD
from newton_raphson import solve_power_flow, calculate_power
from ybus import build_ybus
from run_veragrid import run_veragrid

CASE = "data/case9.m"

off = run_veragrid(CASE, control_q=False)
on = run_veragrid(CASE, control_q=True)

print(f"limits off: converged={off['converged']}, iterations={off['iterations']}")
print(f"limits on:  converged={on['converged']}, iterations={on['iterations']}")
print(f"max |change in V|:     {np.max(np.abs(on['Vm'] - off['Vm'])):.2e} p.u.")
print(f"max |change in angle|: {np.max(np.abs(on['Va_deg'] - off['Va_deg'])):.2e} deg")

# generator Q from my solved voltages: Q injected at the bus + local load Q
case = parse_case(CASE)
mine = solve_power_flow(case)
V = mine["Vm"] * np.exp(1j * np.radians(mine["Va_deg"])) # rebuild complex voltages
S = calculate_power(build_ybus(case), V)
idx = bus_index_map(case["bus"])

print("\n gen bus    Q (MVAr)   Qmin   Qmax   margin to nearest limit")
for g in case["gen"]:
    row = idx[int(g[GEN_BUS])]
    q = S[row].imag * case["baseMVA"] + case["bus"][row, QD] # one generator per bus in case9
    margin = min(g[QMAX] - q, q - g[QMIN])
    print(f"{int(g[GEN_BUS]):8d}   {q:8.2f}   {g[QMIN]:5.0f}  {g[QMAX]:5.0f}   {margin:8.1f}")

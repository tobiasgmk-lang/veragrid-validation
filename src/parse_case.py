"""Read a MATPOWER .m case file into NumPy arrays.

Only reads the three tables a power flow needs (bus, gen, branch) plus baseMVA.
It parses the text; it never runs the file.
"""
import re

import numpy as np

# Column positions (0-based) in MATPOWER's tables, so code can say bus[:, PD]
# instead of bus[:, 2]. Names follow the MATPOWER manual.
BUS_I, BUS_TYPE, PD, QD, GS, BS = 0, 1, 2, 3, 4, 5
VM, VA = 7, 8
PQ, PV, REF = 1, 2, 3                      # bus type codes

GEN_BUS, PG, QG, QMAX, QMIN, VG = 0, 1, 2, 3, 4, 5
GEN_STATUS = 7

F_BUS, T_BUS, BR_R, BR_X, BR_B = 0, 1, 2, 3, 4
TAP, SHIFT, BR_STATUS = 8, 9, 10


def _read_matrix(text, name):
    """Return the numbers between 'mpc.<name> = [' and '];' as a 2-D array."""
    match = re.search(rf"mpc\.{name}\s*=\s*\[(.*?)\];", text, re.S)
    if match is None:
        raise ValueError(f"mpc.{name} not found")
    rows = []
    for line in match.group(1).splitlines():
        line = line.split("%")[0].strip().rstrip(";")   # drop comments and row end
        if line:
            rows.append([float(v) for v in line.split()])
    return np.array(rows)


def parse_case(path):
    with open(path) as f:
        text = f.read()
    base = re.search(r"mpc\.baseMVA\s*=\s*([\d.]+)", text)
    return {
        "baseMVA": float(base.group(1)),
        "bus": _read_matrix(text, "bus"),
        "gen": _read_matrix(text, "gen"),
        "branch": _read_matrix(text, "branch"),
    }

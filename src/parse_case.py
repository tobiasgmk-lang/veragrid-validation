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
    match = re.search(rf"mpc\.{name}\s*=\s*\[(.*?)\];", text, re.S) # re.S: the table can span many lines
    if match is None:
        raise ValueError(f"mpc.{name} not found")
    rows = []
    for line in match.group(1).splitlines():
        line = line.split("%")[0].strip().rstrip(";")   # drop comments and row end
        if line:
            rows.append([float(v) for v in line.split()]) # every number on the row becomes a float
    return np.array(rows) # float array: bus numbers come back as 1.0, 2.0, ... (use int() for lookups)


def bus_index_map(bus):
    """Map each MATPOWER bus number to its row position in the bus table.

    Bus numbers are labels (they can skip, e.g. 1, 2, 5, 10), so code must
    never assume bus k sits in row k-1. Example: {1: 0, 2: 1, 5: 2, 10: 3}.
    """
    return {int(number): row for row, number in enumerate(bus[:, BUS_I])}


def parse_case(path):
    with open(path) as f:
        text = f.read() # read the file as plain text: nothing in it is run
    base = re.search(r"mpc\.baseMVA\s*=\s*([\d.]+)", text) # system base, usually 100 MVA
    return {
        "baseMVA": float(base.group(1)),
        "bus": _read_matrix(text, "bus"),
        "gen": _read_matrix(text, "gen"),
        "branch": _read_matrix(text, "branch"),
    }

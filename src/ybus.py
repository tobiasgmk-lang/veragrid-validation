"""Build the bus admittance matrix (Ybus) from a parsed MATPOWER case.

Each line or transformer is a pi model: series admittance y = 1/(r + jx) between
its two buses, plus half the line charging (j*b/2) to ground at each end.
Taps are not modelled yet (case9 has none), so any off-nominal tap stops the code.
"""
from parse_case import *
import numpy as np

def build_ybus(case):
    bus = case['bus']
    branch = case['branch']
    baseMVA = case['baseMVA']

    n = bus.shape[0] # number of buses, taken from the data

    idx = bus_index_map(bus) # bus number -> row in Ybus

    Ybus = np.zeros((n, n), dtype=complex) # complex, or NumPy silently drops the imaginary part

    for row in branch: # one row per line or transformer
        if row[BR_STATUS] == 0: # skip branches that are out of service
            continue
        assert row[TAP] in (0, 1), f"Off-nominal tap found: {row[TAP]}" # 0 and 1 both mean nominal

        f = idx[int(row[F_BUS])] # row of the from bus
        t = idx[int(row[T_BUS])] # row of the to bus
        r = row[BR_R]
        x = row[BR_X]
        b = row[BR_B] # total line charging, split half to each end

        y = 1.0 / complex(r, x) # series admittance
        b_shunt = 1j * b / 2.0 # line charging at one end

        Ybus[f, f] += y + b_shunt # += not =: several branches meet at one bus
        Ybus[t, t] += y + b_shunt
        Ybus[f, t] += -y # off-diagonals are minus the series admittance
        Ybus[t, f] += -y

    for i in range(n): # fixed bus shunts, given in MW/MVAr at 1 p.u. -> divide by baseMVA
        gs = bus[i, GS]
        bs = bus[i, BS]
        Ybus[i, i] += complex(gs, bs) / baseMVA

    return Ybus

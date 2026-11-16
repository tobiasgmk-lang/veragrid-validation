from parse_case import *
import numpy as np

def build_ybus(case):
    bus = case['bus']
    branch = case['branch']
    baseMVA = case['baseMVA']

    n = bus.shape[0]

    idx = {}
    for i in range(n):
        bus_number = int(bus[i, BUS_I])
        idx[bus_number] = i

    Ybus = np.zeros((n, n), dtype=complex)

    for row in branch:
        if row[BR_STATUS] == 0:
            continue
        assert row[TAP] in (0, 1), f"Non-zero tap ratio found: {row[TAP]}"

        f = idx[int(row[F_BUS])]
        t = idx[int(row[T_BUS])]
        r = row[BR_R]
        x = row[BR_X]
        b = row[BR_B]

        y = 1.0 / complex(r, x)
        b_shunt = 1j * b / 2.0

        Ybus[f, f] += y + b_shunt
        Ybus[t, t] += y + b_shunt
        Ybus[f, t] += -y
        Ybus[t, f] += -y

    for i in range(n):
        gs = bus[i, GS]
        bs = bus[i, BS]
        Ybus[i, i] += complex(gs, bs) / baseMVA
        
    return Ybus
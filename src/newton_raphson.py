import numpy as np

from parse_case import bus_index_map, GEN_BUS, PD, QD, GEN_STATUS, PG

def build_scheduled_power(case):
    bus = case['bus']
    gen = case['gen']
    baseMVA = case['baseMVA']

    n = bus.shape[0]

    idx = bus_index_map(bus)

    S_scheduled = np.zeros(n, dtype=complex) # initialize to complex scheduled array

    S_scheduled -= (bus[:, PD] + 1j * bus[:, QD]) / baseMVA # subtrat the load at buses

    for i in range(gen.shape[0]): # add active generator power
        if gen[i, GEN_STATUS] == 0: # skip generators that are off
            continue
        bus_number = int(gen[i, GEN_BUS])
        bus_index = idx[bus_number]

        S_scheduled[bus_index] += gen[i, PG] / baseMVA 
    return S_scheduled

def calculate_power(Ybus, V):

    S_calc = V * np.conj(Ybus @ V)

    return S_calc
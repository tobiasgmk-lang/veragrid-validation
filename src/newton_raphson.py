"""Newton-Raphson AC power flow, written from scratch to check VeraGrid's.

Unknowns (x): voltage angle at every PV and PQ bus, voltage magnitude at every PQ bus.
Equations (F): P mismatch at every PV and PQ bus, Q mismatch at every PQ bus.
The slack bus has neither: its V and angle are fixed, its P and Q come out at the end.

Sign convention used everywhere in this file:
    F = S_sched - S_calc          (scheduled minus calculated, in p.u.)
    J = d(S_calc) / dx            (Jacobian of the CALCULATED power)
    dx = solve(J, F)              (so x_new = x + dx moves S_calc towards S_sched)
"""
import numpy as np

from parse_case import (bus_index_map, BUS_I, BUS_TYPE, PD, QD, PV, PQ, REF,
                        GEN_BUS, GEN_STATUS, PG, VG)
from ybus import build_ybus


def build_scheduled_power(case):
    bus = case['bus']
    gen = case['gen']
    baseMVA = case['baseMVA']

    n = bus.shape[0]

    idx = bus_index_map(bus)

    S_scheduled = np.zeros(n, dtype=complex) # initialize to complex scheduled array

    S_scheduled -= (bus[:, PD] + 1j * bus[:, QD]) / baseMVA # subtract the load at buses

    for i in range(gen.shape[0]): # add active generator power
        if gen[i, GEN_STATUS] == 0: # skip generators that are off
            continue
        bus_number = int(gen[i, GEN_BUS])
        bus_index = idx[bus_number]

        S_scheduled[bus_index] += gen[i, PG] / baseMVA # Q left out: unknown at PV buses
    return S_scheduled


def calculate_power(Ybus, V):

    S_calc = V * np.conj(Ybus @ V) # power injected at each bus for voltages V

    return S_calc


def bus_types(case):
    bus = case['bus']

    pv = np.where(bus[:, BUS_TYPE] == PV)[0] # row positions of PV buses
    pq = np.where(bus[:, BUS_TYPE] == PQ)[0] # row positions of PQ buses
    slack = np.where(bus[:, BUS_TYPE] == REF)[0]
    assert len(slack) == 1, f"expected one slack bus, found {len(slack)}"

    pvpq = np.concatenate((pv, pq)) # every bus with an unknown angle

    return pv, pq, slack, pvpq


def initial_voltage(case):
    bus = case['bus']
    gen = case['gen']
    idx = bus_index_map(bus)

    Vm = np.ones(bus.shape[0]) # PQ buses start at 1 p.u.
    for i in range(gen.shape[0]): # PV and slack buses start (and stay) at the setpoint Vg
        if gen[i, GEN_STATUS] == 0:
            continue
        Vm[idx[int(gen[i, GEN_BUS])]] = gen[i, VG]

    return Vm.astype(complex) # all angles 0: V = Vm at 0 degrees


def mismatch(S_sched, S_calc, pq, pvpq):

    delta_S = S_sched - S_calc

    F = np.concatenate((
        delta_S[pvpq].real, # P mismatch at PV and PQ buses
        delta_S[pq].imag    # Q mismatch at PQ buses only
    ))
    return F


def build_jacobian(Ybus, V, pq, pvpq):
    """YOURS TO WRITE. Return the (len(pvpq) + len(pq)) square Jacobian of S_calc.

    Rows follow F:      [P at pvpq ; Q at pq]
    Columns follow dx:  [angle at pvpq ; |V| at pq]
    so J = [[dP/dtheta, dP/d|V|],
            [dQ/dtheta, dQ/d|V|]]
    """
    raise NotImplementedError("build_jacobian")


def update_voltage(V, dx, pq, pvpq):
    """YOURS TO WRITE. Apply the Newton step dx and return the new complex V.

    dx is ordered like the Jacobian columns: [angle at pvpq ; |V| at pq].
    """
    raise NotImplementedError("update_voltage")


def solve_power_flow(case, tol=1e-6, max_iter=25):
    Ybus = build_ybus(case)
    S_sched = build_scheduled_power(case)
    pv, pq, slack, pvpq = bus_types(case)
    V = initial_voltage(case)

    converged = False
    iterations = 0 # counts Newton steps taken (Jacobian solves)
    while True:
        S_calc = calculate_power(Ybus, V)
        F = mismatch(S_sched, S_calc, pq, pvpq)
        error = np.max(np.abs(F)) # largest mismatch in p.u., same test VeraGrid uses

        if error < tol: # converged
            converged = True
            break
        if iterations == max_iter: # give up, report failure
            break

        J = build_jacobian(Ybus, V, pq, pvpq)
        dx = np.linalg.solve(J, F) # Newton step: J dx = F
        V = update_voltage(V, dx, pq, pvpq)
        iterations += 1

    S_calc = calculate_power(Ybus, V)
    return {
        "bus": [int(b) for b in case['bus'][:, BUS_I]],
        "Vm": np.abs(V),
        "Va_deg": np.degrees(np.angle(V)),
        "converged": converged,
        "iterations": iterations,
        "error": error,
        # sum of every bus injection = what the network consumed = losses
        "losses_MW": float(np.sum(S_calc.real) * case['baseMVA']),
    }

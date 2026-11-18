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
    # calculate current and unit voltage phasor
    I = Ybus @ V # current injected at each bus
    E = V / np.abs(V) # V scaled to magnitude 1: keeps only the direction (angle)

    # calculate the partial derivative of complex power
    # (matrix form of the product rule on S = V * conj(Ybus V), from MATPOWER's technical note 2)
    dS_dtheta = ( # how S at every bus changes when each angle changes
    1j * np.diag(V)
    @ np.conj(np.diag(I) - Ybus @ np.diag(V))
    )

    dS_dVm = ( # how S at every bus changes when each |V| changes
    np.diag(V) @ np.conj(Ybus @ np.diag(E))
    + np.conj(np.diag(I)) @ np.diag(E)
    )

    # 4 jacobian submatrices: real part = P rows, imaginary part = Q rows
    # (angles and |V| are real numbers, so d(Re S) = Re(dS))
    H = dS_dtheta[np.ix_(pvpq, pvpq)].real # dP/dtheta
    N = dS_dVm[np.ix_(pvpq, pq)].real      # dP/d|V|  (only pq buses have an unknown |V|)
    M = dS_dtheta[np.ix_(pq, pvpq)].imag   # dQ/dtheta (only pq buses have a Q equation)
    L = dS_dVm[np.ix_(pq, pq)].imag        # dQ/d|V|

    # assemble the full jacobian: rows follow F, columns follow dx. the slack is in neither
    J = np.block([
        [H, N],
        [M, L]
    ])

    return J


def update_voltage(V, dx, pq, pvpq):
    # split the newton step into angle and magnitude components
    # (same order as the jacobian columns: angles first, then |V|)
    n_angle = len(pvpq)
    dtheta = dx[:n_angle] # radians
    dVm = dx[n_angle:] # p.u.

    # extract current voltage magnitude and angle
    Vm = np.abs(V)
    theta = np.angle(V)

    # apply corrections to the voltage magnitude and angle
    Vm[pq] += dVm # only pq buses: pv and slack |V| stay at their setpoint Vg
    theta[pvpq] += dtheta # every bus except the slack, which stays at 0

    # reconstruct the complex voltage phasor: V = |V| e^(j theta)
    V = Vm * np.exp(1j * theta)

    return V 


def solve_power_flow(case, tol=1e-6, max_iter=25):
    Ybus = build_ybus(case)
    S_sched = build_scheduled_power(case)
    pv, pq, slack, pvpq = bus_types(case)
    V = initial_voltage(case)

    converged = False
    iterations = 0 # counts Newton steps taken (Jacobian solves)
    history = [] # largest mismatch before each step, to show how fast NR converges
    while True:
        S_calc = calculate_power(Ybus, V)
        F = mismatch(S_sched, S_calc, pq, pvpq)
        error = np.max(np.abs(F)) # largest mismatch in p.u., same test VeraGrid uses
        history.append(error)

        if error < tol: # converged. smaller tol = more accurate answer, maybe one more step
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
        "mismatch_history": history,
        # sum of every bus injection = what the network consumed = losses
        "losses_MW": float(np.sum(S_calc.real) * case['baseMVA']),
    }

# VeraGrid power flow validation (IEEE 9-bus)

## What this is and why

This is me trying to build a power flow so that I can learn the mechanics behind it. It's about being able to trust a tool that I would test, and to identify mistakes in it as well.

I did this because I've been fascinated by this topic. I want an internship in it and I would like to work in it, but I have no experience in it.

The project validates VeraGrid's AC power flow against a Newton–Raphson solver that I wrote myself, and tests and documents the comparison. Both solvers read the same IEEE 9-bus case file (`case9.m` from MATPOWER), so any difference comes from the method or a setting, not the data. I also ran MATPOWER 8.0 (in MATLAB) as a third check.

## How it works

### What a power flow does

A power flow calculates the power from each generator, the power consumed by each load and the power flowing in each line. It finds P, Q, V and θ at every bus when only some of them are known, usually the scheduled power. The system is nonlinear: changing the power at one bus changes the voltage at every other bus, so it has to be solved iteratively.

Each bus has two known and two unknown quantities:

| Bus type | Known | Unknown |
|---|---|---|
| PQ (load) | P, Q | \|V\|, θ |
| PV (generator) | P, \|V\| | Q, θ |
| Slack (reference) | \|V\|, θ | P, Q |

Once the voltages are known, everything else follows: power injections at every bus, flows in every line and the losses.

### The case file

`data/case9.m` is plain text with four parts:

- `baseMVA`: the system base (100 MVA). Powers in the file are in MW/MVAr and are divided by this to get per-unit.
- `bus`: one row per bus. Bus number, type (1 = PQ, 2 = PV, 3 = slack), load (Pd, Qd) and any shunt to ground (Gs, Bs).
- `gen`: one row per generator. Which bus it's on, its scheduled P (Pg), its voltage setpoint (Vg) and its reactive limits.
- `branch`: one row per line or transformer. From and to bus, r, x, the line charging b, and the tap ratio.

A note on line charging: a long line is almost like a long, thin capacitor. It generates reactive power just by being energised, which raises the voltage. In case9 the lines produce about 140 MVAr, more than the 115 MVAr the loads consume, which is why most buses sit above 1.0 p.u. The π model is where that capacitance sits: half of `b` at each end of the line.

### The steps

1. **Build the electrical network** (`ybus.py`). Each branch adds its series admittance and half its line charging to Ybus. I check this matrix entry by entry against VeraGrid's before running any power flow.
2. **Establish the scheduled power injections** (`build_scheduled_power`): generation minus load at each bus.
3. **Identify what is known and unknown** (`bus_types`).
4. **Make an initial estimate of the voltages** (`initial_voltage`): every angle at 0, PQ buses at 1 p.u., generator buses at their setpoint.
5. **Calculate the power from those voltages** (`calculate_power`): Ybus gives the currents, and the currents give the power injections.
6. **Compare the calculated and scheduled power** (`mismatch`).
7. **Find out how sensitive the power is to each voltage** (`build_jacobian`).
8. **Correct the voltages** (`np.linalg.solve`, then `update_voltage`), and repeat from step 5 until the largest mismatch is below 1e-6 p.u.
9. **Compare and test the results** against VeraGrid (`tests/`).

The Jacobian has four blocks:

| Block | Meaning |
|---|---|
| H | how active power changes when voltage angles change |
| N | how active power changes when voltage magnitudes change |
| M | how reactive power changes when voltage angles change |
| L | how reactive power changes when voltage magnitudes change |

The error drops very quickly near the answer because of quadratic convergence. Each step roughly squares the previous error:

| Step | Largest mismatch (p.u.) |
|---|---|
| 0 | 1.63 |
| 1 | 0.19 |
| 2 | 0.0021 |
| 3 | 0.00000034 |

## How to run it

You need Python 3.13 and git.

```
git clone https://github.com/tobiasgmk-lang/veragrid-validation.git
cd veragrid-validation
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Then, from the project folder:

| Command | What it does |
|---|---|
| `python -m pytest -v -s` | runs both tests and writes `results/case9_comparison.csv` and `results/case9_summary.csv` |
| `python src/run_veragrid.py data/case9.m` | prints VeraGrid's voltages by bus |
| `python scripts/check_q_limits.py` | runs VeraGrid with reactive limits off and on, and shows each generator's margin |

Optional, MATPOWER (in MATLAB, from the project folder):

```matlab
mpopt = mpoption('pf.tol', 1e-6);
r = runpf('case9', mpopt);
writematrix(r.bus(:, [1 8 9]), 'results/matpower_case9.csv');
```

`data/case14.m` is included but not supported yet: it has transformer taps, and the code stops on purpose when it finds one.

The case files in `data/` are from [MATPOWER](https://github.com/MATPOWER/matpower) and are used under MATPOWER's BSD 3-clause license.

## Results

Case9, after matching buses by bus number. Differences are the largest over all 9 buses.

| | Mine | VeraGrid 6.8.4 | MATPOWER 8.0 |
|---|---|---|---|
| Converged | yes | yes | yes |
| Iterations | 3 | 3 | 3 |
| Total losses (MW) | 4.641023 | 4.641023 | 4.641023 |
| Largest \|ΔV\| vs mine (p.u.) | – | 4.4e-16 | 3.3e-15 |
| Largest Δθ vs mine (degrees) | – | 8.9e-15 | 1.1e-14 |

The Ybus matrices (all 81 entries) differ by at most 4.4e-16 p.u.

## What differed and why

Nothing differed apart from floating-point rounding. The voltage magnitudes over the 9 buses agree with VeraGrid to 4.4e-16 p.u., the angles to 8.9e-15°, and the total losses to 2.9e-14 MW. MATPOWER agrees with mine to 3.3e-15 p.u. and 1.1e-14°. All three solvers took 3 iterations.

Things that could have caused differences, and what I did about them:

- **Bus order.** I match the results by bus number and rearrange them, rather than relying on array position.
- **Reactive limits.** These are off by default in VeraGrid. I turned them on to check, and nothing changed. The limits are ±300 MVAr, which leaves a lot of headroom: the closest generator was 273.0 MVAr from its limit.

Bugs the tests caught in my own code:

- **An indentation bug in Ybus.** Only the last branch was being added, so Ybus had 4 non-zero entries instead of 27.
- **The bus number used as an array position.** The code crashed at bus 9. If there had been no bus 9, it wouldn't have crashed at all, and every bus would have been silently compared with the wrong one.

## Limitations

- Three solvers agreeing to the last digit is reassuring, but it's also suspicious. They have the same formulation, so they can't catch a mistake that they'd all make.
- It's only a small case. A 9-bus system converges easily.
- It doesn't test transformer taps, bus shunts, out-of-service branches, reactive limits, or several generators at one bus, because case9 has none of these.
- The tolerances assume that the smallest real modelling error is about 1e-3 p.u. That's an estimate, not a measurement.

## Next steps

- **case14:** add taps to the Ybus code and test the bus shunt (bus 9).
- **Break the model on purpose**, for example by removing b/2 from one line, to see the test fail and to measure the real error the tolerances rely on.
- **Add reactive limits to my solver.**

## How I used AI

During the creation of this project I used Claude (Opus 5.5, Anthropic), and another AI assistant and Stack Overflow for parts of the Ybus, the Jacobian and the test comparison. Whilst I would like to be at the stage where I could do this myself, I tried to understand, the best I could, what the AI was doing.

I wrote:
- the Ybus
- the scheduled power
- the bus types
- the mismatch vector
- the Jacobian and voltage updates (the matrix formulas for the Jacobian came from Claude, from MATPOWER's technical note)
- the comparison logic for both tests

The AI wrote:
- the parser
- the VeraGrid runner
- the CSV writer
- the reactive (Q) limit script
- the solver loop and the starting voltages
- most of the comments
- this README, assembled from what I dictated

In order to build my understanding I tried to predict results before running them, explain the structure back as I understood it, check independently with VeraGrid and MATPOWER, and test each part as it was built.

I've learnt Python through coursework. I would not call myself an expert. I took notes the best I could on my understanding throughout this, and with my notes I could rebuild the `build_jacobian` function. I will try to do it later. I also used the AI to make sure that my code only considered things that were relevant and was modular, so that I could use it in other cases. I tried to make it as efficient as possible and easy for other people to follow along with what is happening.

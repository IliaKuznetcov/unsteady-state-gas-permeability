# FDM Verification Notes

## Scope

These notes summarize the numerical verification of the explicit finite-difference
implementation of the ideal-gas Darcy–Klinkenberg pulse-decay model.

The mathematical formulation and required verification checks are documented in:

[Finite-Difference Formulation of a Gas Pulse-Decay Experiment](Equation_derivations_FDM_USS.pdf)

The verified reference implementation uses the **storage-consistent boundary
formulation** and the transformed pressure variable

```math
\phi = (P+b)^2
```

with explicit Forward Euler time integration.

The purpose of these notes is to record what has actually been tested in the
Python implementation and to distinguish numerical verification from later
comparison with experimental data.

---

## 1. Baseline Reference Case

The current baseline problem uses:

| Parameter | Value |
|---|---:|
| Upstream volume, $V_u$ | $5.0\times10^{-6}\ \mathrm{m^3}$ |
| Downstream volume, $V_d$ | $5.0\times10^{-6}\ \mathrm{m^3}$ |
| Core length, $L$ | $0.05\ \mathrm{m}$ |
| Core diameter, $D$ | $0.0254\ \mathrm{m}$ |
| Porosity, $\varepsilon$ | 0.10 |
| Intrinsic permeability, $k_\ell$ | $1.0\times10^{-17}\ \mathrm{m^2}$ |
| Gas viscosity, $\mu$ | $1.8\times10^{-5}\ \mathrm{Pa\,s}$ |
| Klinkenberg factor, $b$ | $1.0\times10^5\ \mathrm{Pa}$ |
| Initial upstream pressure | $1.2\ \mathrm{MPa}$ |
| Initial downstream/core pressure | $1.0\ \mathrm{MPa}$ |
| Spatial intervals, $N$ | 40 |
| Simulation duration | $5000\ \mathrm{s}$ |
| Stability safety factor | 0.8 |

For $N=40$,

```math
\Delta x = 1.25\times10^{-3}\ \mathrm{m}
```

and the initial stability-based timestep is approximately

```math
\Delta t = 8.653846\times10^{-2}\ \mathrm{s}.
```

The baseline calculation gives approximately:

| Quantity | Result |
|---|---:|
| Final upstream pressure | $1.080270\ \mathrm{MPa}$ |
| Final downstream pressure | $1.080269\ \mathrm{MPa}$ |
| Final differential pressure | $0.001\ \mathrm{kPa}$ |
| Analytical continuous equilibrium | $1.079786\ \mathrm{MPa}$ |
| Discrete initial-inventory equilibrium | $1.080291\ \mathrm{MPa}$ |
| Error vs. analytical equilibrium | $+0.484\ \mathrm{kPa}$ |
| Error vs. discrete equilibrium | $-0.022\ \mathrm{kPa}$ |
| Maximum relative inventory drift | $2.012071\times10^{-5}$ |

The tank-pressure histories show the expected qualitative behavior:

```math
P_u(t)\downarrow P_{\mathrm{eq}},
\qquad
P_d(t)\uparrow P_{\mathrm{eq}},
\qquad
\Delta P(t)\rightarrow 0.
```

---

## 2. Inventory Diagnostic

The discrete pressure-volume inventory is evaluated on the node-centered grid
using the trapezoidal rule:

```math
I^n
=
V_u P_0^n
+
V_d P_N^n
+
\varepsilon S\Delta x
\left[
\frac{1}{2}P_0^n
+
\sum_{i=1}^{N-1}P_i^n
+
\frac{1}{2}P_N^n
\right].
```

The relative inventory drift is

```math
E_I^n = \frac{I^n-I^0}{I^0}.
```

This quantity is a diagnostic rather than an exactly conserved discrete
invariant for the transformed FDM. The physical storage is linear in pressure,
whereas the numerical update advances $\phi=(P+b)^2$.

For the $N=40$ baseline calculation, the maximum relative drift is about

```math
2.01\times10^{-5},
```

or roughly 20 ppm.

---

## 3. Initial-Condition Discretization Bias

The idealized continuous initial condition assigns the higher upstream pressure
to the tank/core interface at $x=0$, while the core for $x>0$ begins at the
lower downstream pressure.

On the node-centered grid, the trapezoidal inventory assigns half of the first
core interval's pore storage to node 0. The discrete initial state therefore
contains slightly more gas than the corresponding continuous initial condition.

The leading inventory bias is approximately

```math
\Delta I_{\mathrm{grid}}
\approx
\frac{\varepsilon S\Delta x}{2}
\left(P_{u,i}-P_{d,i}\right),
```

so the associated equilibrium-pressure error is $O(\Delta x)$.

For the $N=40$ baseline case, the discrete equilibrium is about
$0.505\ \mathrm{kPa}$ above the continuous analytical equilibrium. The small
negative inventory drift during the transient offsets about $0.022\ \mathrm{kPa}$
of this bias, leaving the numerical equilibrium about $0.484\ \mathrm{kPa}$
above the continuous analytical value.

This is primarily an **initial-condition/grid representation error**, not
time-accumulated inventory drift.

---

## 4. Grid Refinement

The grid-refinement study uses

```math
N = 10,\ 20,\ 40,\ 80,
```

with the timestep reduced consistently with the explicit stability condition.

### Results

| $N$ | $\Delta x$ [m] | $\Delta t$ [s] | Error vs. analytical equilibrium [kPa] | Error vs. discrete equilibrium [kPa] | Max. inventory error |
|---:|---:|---:|---:|---:|---:|
| 10 | $5.000000\times10^{-3}$ | 1.384615 | 1.934667 | -0.086739 | $8.017994\times10^{-5}$ |
| 20 | $2.500000\times10^{-3}$ | 0.346154 | 0.967295 | -0.043408 | $4.016287\times10^{-5}$ |
| 40 | $1.250000\times10^{-3}$ | 0.0865385 | 0.483615 | -0.021736 | $2.012071\times10^{-5}$ |
| 80 | $6.250000\times10^{-4}$ | 0.0216346 | 0.241793 | -0.010882 | $1.007600\times10^{-5}$ |

For an asymptotic error model

```math
E(\Delta x)\approx C\,\Delta x^p,
```

the observed order is estimated from two successive meshes as

```math
p
=
\frac{
\ln\left(|E_{\mathrm{coarse}}|/|E_{\mathrm{fine}}|\right)
}{
\ln\left(\Delta x_{\mathrm{coarse}}/\Delta x_{\mathrm{fine}}\right)
}.
```

The observed orders are approximately:

| Refinement | $p$, analytical equilibrium error | $p$, discrete equilibrium error | $p$, inventory drift |
|---|---:|---:|---:|
| 10 → 20 | 1.0001 | 0.9987 | 0.9974 |
| 20 → 40 | 1.0001 | 0.9979 | 0.9972 |
| 40 → 80 | 1.0001 | 0.9981 | 0.9978 |

The coupled scheme therefore exhibits approximately **first-order global
convergence** for these metrics.

This does not contradict the second-order centered interior approximation.
The observed global behavior includes the endpoint treatment, Forward Euler
time integration, the discontinuous initial condition, and the $O(\Delta x)$
initial storage representation.

---

## 5. Time Refinement

The temporal-refinement study fixes $N=40$ and compares successive timestep
factors

```math
1,\quad
\frac{1}{2},\quad
\frac{1}{4},\quad
\frac{1}{8},\quad
\frac{1}{16}.
```

All solutions are compared at the same physical time,

```math
t^* = 500.019231\ \mathrm{s},
```

by choosing integer step counts that align exactly across refinement levels.

### Transient results

| Timestep factor | $\Delta t$ [s] | $P_u$ [MPa] | $P_d$ [MPa] | $\Delta P$ [kPa] |
|---:|---:|---:|---:|---:|
| 1.0000 | $8.653846\times10^{-2}$ | 1.107409 | 1.053032 | 54.376813 |
| 0.5000 | $4.326923\times10^{-2}$ | 1.107426 | 1.053039 | 54.387207 |
| 0.2500 | $2.163462\times10^{-2}$ | 1.107433 | 1.053042 | 54.391691 |
| 0.1250 | $1.081731\times10^{-2}$ | 1.107437 | 1.053043 | 54.393818 |
| 0.0625 | $5.408654\times10^{-3}$ | 1.107438 | 1.053044 | 54.394856 |

Because no exact transient solution is available, temporal order is estimated
by self-convergence. For three successive solutions,

```math
D_1 = |Q_{\Delta t}-Q_{\Delta t/2}|,
\qquad
D_2 = |Q_{\Delta t/2}-Q_{\Delta t/4}|,
```

and

```math
p_t
=
\frac{\ln(D_1/D_2)}{\ln 2}.
```

### Observed temporal order

| Refinement triplet | $p_t(P_u)$ | $p_t(P_d)$ | $p_t(\Delta P)$ |
|---|---:|---:|---:|
| 1 → 1/2 → 1/4 | 1.2880 | 1.4055 | 1.2128 |
| 1/2 → 1/4 → 1/8 | 1.1071 | 1.1598 | 1.0764 |
| 1/4 → 1/8 → 1/16 | 1.0479 | 1.0734 | 1.0337 |

The measured temporal order approaches

```math
p_t = 1,
```

consistent with the expected asymptotic first-order accuracy of Forward Euler.

---

## 6. Automated Verification Tests

The reusable solver in `fdm_solver.py` is exercised directly by
`fdm_verification_tests.py`.

The current automated tests are:

### Uniform-pressure test

Set

```math
P_{u,i}=P_{d,i}.
```

With no pressure gradient, every node and both tank pressures remain unchanged
to numerical roundoff.

**Status: PASS**

### Zero-permeability test

Set

```math
k_\ell=0.
```

The normal stability timestep expression is bypassed because its limiting value
is infinite. No pressure state changes.

**Status: PASS**

### Baseline regression test

The reusable `run_fdm()` implementation is required to reproduce the previously
verified $N=40$, 5000 s baseline solution after refactoring.

**Status: PASS**

---

## 7. Verification Status

The current storage-consistent reference implementation has demonstrated:

- correct qualitative pressure evolution;
- positivity and explicit coefficient checks;
- approach to equilibrium;
- pressure-volume inventory diagnostics;
- grid convergence;
- first-order global convergence for the tested equilibrium/inventory metrics;
- temporal self-convergence toward first-order Forward Euler behavior;
- uniform-state preservation;
- zero-permeability behavior;
- regression consistency after extraction into a reusable solver.

The following checks listed in the original technical notes have **not yet been
used as acceptance criteria for the current reference implementation**:

- direct-boundary versus storage-consistent boundary comparison;
- deliberately unstable timestep demonstration.

These are useful diagnostic studies, but the verified forward model used for
subsequent parameter estimation is the storage-consistent formulation.

---

## 8. Interpretation of Accuracy Claims

Several numerical properties should remain conceptually separate:

- **stability** concerns whether perturbations remain bounded under the explicit
  timestep restrictions;
- **monotonicity/convexity** concerns whether one-step updates introduce new
  extrema or overshoot;
- **inventory conservation** is assessed independently through $E_I(t)$;
- **verification** asks whether the numerical implementation solves the stated
  equations consistently;
- **validation** asks whether the mathematical model reproduces experimental
  pulse-decay measurements.

The centered interior second derivative is formally second order in space, but
the order of the complete coupled calculation must be established by refinement
of the full model rather than inferred from that stencil alone.

---

## 9. Next Stage

The next development stage is experimental parameter estimation.

The verified FDM solver will be used as the forward model inside a nonlinear
least-squares problem comparing simulated and measured upstream/downstream
pressure histories.

The initial inverse-model parameters are expected to be

```math
\theta = (k_\ell,\ b),
```

while calibrated tank volumes, dead volumes, core dimensions, viscosity, and
porosity are treated as measured model inputs where possible.

After the FDM fitting workflow is established on real pulse-decay data, a
cell-centered finite-volume formulation with shared face fluxes and
pressure-based storage can be developed and compared against the FDM reference.

# Unsteady-State Gas Permeability

Numerical modeling and parameter-estimation tools for unsteady-state gas
permeability measurements using pulse-decay experiments.

The project currently contains a verified explicit finite-difference model
for one-dimensional ideal-gas flow through a porous core with
Darcy–Klinkenberg permeability.

The governing equations and numerical formulation are derived in:

`Equation_derivations_FDM_USS.pdf`

---

## Physical Model

The model describes gas flow between upstream and downstream storage volumes
connected through a cylindrical porous core.

Gas density is assumed to follow the ideal-gas relation

```math
\rho = C_g P
```

and apparent gas permeability is represented using the Klinkenberg relation

```math
k_a(P) = k_\ell \left(1 + \frac{b}{P}\right)
```

where:

- $k_\ell$ is the intrinsic/liquid-equivalent permeability;
- $b$ is the Klinkenberg slip factor;
- $P$ is gas pressure.

Combining mass conservation with Darcy flow gives

```math
\varepsilon
\frac{\partial P}{\partial t}
=
\frac{k_\ell}{\mu}
\frac{\partial}{\partial x}
\left[
(P+b)
\frac{\partial P}{\partial x}
\right]
```

The pressure transformation

```math
\phi = (P+b)^2
```

is used to obtain the explicit finite-difference formulation implemented in
this repository.

---

## Current FDM Implementation

The current solver uses:

- a one-dimensional node-centered spatial grid;
- explicit Forward Euler time integration;
- centered second-order spatial differences for interior nodes;
- pressure-dependent Darcy–Klinkenberg transport coefficients;
- storage-consistent upstream and downstream tank boundary conditions;
- an automatically calculated explicit stability timestep;
- runtime stability, positivity, and finite-value checks;
- pressure-volume inventory diagnostics.

The reusable forward solver is implemented in `fdm_solver.py` and can be
called from other scripts using:

```python
from fdm_solver import run_fdm
```

---

## Numerical Verification

The FDM implementation has been checked using independent limiting-case,
regression, grid-refinement, and time-refinement tests.

### Limiting cases

The automated verification suite confirms that:

- a uniform initial pressure remains unchanged;
- zero permeability produces no pressure evolution;
- the reusable solver reproduces the previously verified baseline solution.

Run the verification suite with:

```bash
python fdm_verification_tests.py
```

### Grid refinement

Spatial refinement has been tested using

```math
N = 10,\ 20,\ 40,\ 80
```

The coupled FDM solution shows approximately first-order global convergence
for equilibrium-pressure and inventory-error metrics.

This first-order behavior describes the complete coupled formulation,
including the tank boundary treatment and initial-condition representation,
rather than only the centered interior finite-difference stencil.

Run with:

```bash
python fdm_grid_refinement.py
```

### Time refinement

Temporal refinement has been tested using successive timestep reductions

```math
\Delta t,\quad
\frac{\Delta t}{2},\quad
\frac{\Delta t}{4},\quad
\frac{\Delta t}{8},\quad
\frac{\Delta t}{16}
```

The observed temporal convergence order approaches

```math
p_t \approx 1
```

consistent with the expected first-order accuracy of Forward Euler time
integration.

Run with:

```bash
python fdm_time_refinement.py
```

---

## Repository Structure

```text
fdm_solver.py
    Reusable explicit FDM forward solver.

fdm_baseline.py
    Baseline pulse-decay simulation, plotting, equilibrium diagnostics,
    and pressure-volume inventory analysis.

fdm_grid_refinement.py
    Spatial grid-refinement and convergence study.

fdm_time_refinement.py
    Temporal refinement and self-convergence study.

fdm_verification_tests.py
    Uniform-pressure, zero-permeability, and baseline regression tests.

Equation_derivations_FDM_USS.pdf
    Mathematical derivation and technical notes for the FDM formulation.
```

---

## Verified FDM Milestone

The verified forward-model implementation is tagged as:

`fdm-v0.1-verified`

This tag provides a fixed reference point before introducing experimental
data handling and inverse parameter estimation.

---

## Current Limitations

The current FDM formulation is intended primarily as a verified reference
model.

Because pressure is represented on a node-centered grid, the initial pressure
discontinuity between the upstream reservoir and the core introduces a
grid-dependent pressure-volume inventory bias of order $O(\Delta x)$.

The transformed-variable explicit FDM also does not conserve the discrete
pressure-volume inventory exactly, although the inventory drift decreases
under grid refinement.

These behaviors are quantified by the verification scripts rather than
hidden from the analysis.

---

## Planned Development

The next stages of the project are:

1. import calibrated upstream/downstream tank and dead-volume data;
2. load experimental upstream and downstream pressure histories;
3. use the verified FDM model as a forward model in a nonlinear
   least-squares inverse problem;
4. estimate permeability $k_\ell$ and Klinkenberg factor $b$ from
   pulse-decay experiments;
5. compare modeled and measured pressure histories and residuals;
6. develop a conservative finite-volume formulation using shared face fluxes;
7. compare FDM and FVM parameter estimates and conservation behavior;
8. build a PySide6 interface for experimental analysis and visualization.

---

## Purpose

The long-term objective is to develop a transparent and reproducible
unsteady-state gas-permeability analysis workflow in which the governing
equations, numerical discretization, verification studies, experimental data
processing, and parameter estimation remain explicitly connected.

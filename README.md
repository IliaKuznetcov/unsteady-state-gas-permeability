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

\[
\rho = C_g P
\]

and apparent gas permeability is represented using the Klinkenberg relation

\[
k_a(P) = k_\ell \left(1 + \frac{b}{P}\right),
\]

where:

- \(k_\ell\) is the intrinsic/liquid-equivalent permeability;
- \(b\) is the Klinkenberg slip factor;
- \(P\) is gas pressure.

Combining mass conservation with Darcy flow gives

\[
\varepsilon
\frac{\partial P}{\partial t}
=
\frac{k_\ell}{\mu}
\frac{\partial}{\partial x}
\left[
(P+b)
\frac{\partial P}{\partial x}
\right].
\]

The transformation

\[
\phi = (P+b)^2
\]

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

The reusable forward solver is implemented in:

`fdm_solver.py`

and can be called from other scripts using:

```python
from fdm_solver import run_fdm
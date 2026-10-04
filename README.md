# Unsteady-State Gas Permeability — FDM

Finite-difference implementation of the ideal-gas Darcy–Klinkenberg
pulse-decay model.

The implementation follows the derivation in:

`Equation_derivations_FDM_USS.pdf`

Initial objectives:

- implement and verify the explicit FDM formulation;
- compare direct and storage-consistent boundary schemes;
- add numerical verification tests;
- later develop a PySide6 interface.
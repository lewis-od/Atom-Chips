# Atom Chips
Numerical computations of the magnetic trap profiles of various current-carrying wire configurations, used for trapping ultracold atoms on a chip.

[U Wire](U%20Wire/) and [Z Wire](Z%20Wire/) use the Biot-Savart law and the finite-element method to calculate the magnetic field for a finite, 1d, current carrying wire in the U and Z configurations.

[Surface](Surface/) uses the MATLAB PDE toolkit to calculate the current density on a 2D conductor with a given geometry, then calculates the magnetic field associated with this at a given distance from the surface.

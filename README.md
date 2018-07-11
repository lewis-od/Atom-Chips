# Atom Chips
Numerical computations of the magnetic trap profiles of various current-carrying wire configurations, used for trapping ultracold atoms on a chip.

[U Wire](U%20Wire/) and [Z Wire](Z%20Wire/) use the Biot-Savart law and the finite-element method to calculate the magnetic field for a finite, 1d, current carrying wire in the U and Z configurations.

[Surface 2.0](Surface%202.0/) uses the MATLAB PDE toolkit to calculate the current density on a conducting channel in a two-dimensional electron gas with a given geometry (in this case, a 'Z' shape), then calculates the magnetic field associated with this. The screened electrostatic potential due to ionised donors can optionally be included. Properties of the magentic trap formed, such as trap depth, frequency, and loss rate due to spin-flips are then calculated. See [this paper](https://arxiv.org/abs/1708.01184) for more information about this type of atom chip. Appendix A explains the calculation of the screened electrostatic potential.

[Surface](Surface/) contains code for a (failed) attempt to calculate the magnetic field due to a 2D current distribution using Fourier transforms

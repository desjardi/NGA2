# Liquid slab in a Couette flow

- This case simulates the deformation of a liquid slab in a Couette flow to better understand how to:
    - properly implement sliding walls and Navier slip boundary conditions
    - assess how to model contact line physics reliably with semi-Lagrangian VOF
- The case assumes:
    - The wall velocity is $$U_w=1$$
    - The channel height is $$h=1$$
- The specific quantities to set in the `input` file are:
    - `Re`, the liquid Reynolds number, which is expressed as $$\mathrm{Re} = \dfrac{h U_w}{\mu_l}$$
    - `Oh`, the liquid Ohnesorge number, which is given by $$\mathrm{Oh} = \dfrac{\mu_l}{\sqrt{\rho_l h \sigma}}$$
    - `R`, the density ratio, defined by $$R=\dfrac{\rho_l}{\rho_g}$$
    - and `M`, the viscosity ratio, given by $$M=\dfrac{\mu_l}{\mu_g}$$
- Finally, the contact angle is specified in degrees using the `input` file variable `CA`

<p align="center"><img src="./slab.png"/></p>
<p align="center">Snapshot of the simulation.</p>

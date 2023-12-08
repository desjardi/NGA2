# Liquid film retraction in gas

- This is a two-phase flow simulation of an unstable retracting liquid film in (potentially) turbulent air, following some of the work in [Mirjalili S., Chan W. H. R., Mani A. (2018) High Fidelity Simulations of Micro-Bubble Shedding from Retracting Thin Gas Films in the Context of Liquid-Liquid Impact, 32nd Symposium on Naval Hydrodynamics](https://arxiv.org/pdf/1811.12352)
- The mesh is set in the `input` file based on:
    - The number of cells, given by `nx`, `ny`, and `nz`.
    - The domain size adimensionalized by the film thickness, given by `Lx`, `Ly`, and `Lz`.
    - Note that at this point, the mesh is uniform Cartesian.
- The main adimensional physical parameters to set in the `input` file are:
    - The `Ohnesorge number` defined as
    $$\mathrm{Oh}=\dfrac{\mu_l}{\sqrt{\rho_l h \sigma}}$$.
    - The `Density ratio`, defined as
    $$r=\dfrac{\rho_l}{\rho_g}$$.
    - The `Viscosity ratio`, defined as
    $$m=\dfrac{\mu_l}{\mu_g}$$.
- The `Inflow velocity` can also *(optionally)* be specified, by default it is set to the adimensional Taylor--Culick retraction velocity of $\sqrt{2}$.

<p align="center"><img src="./snap_t2.png"/></p>
<p align="center">Snaphost of the simulation at $t = 2\tau_L$ (after two eddy-turnover times).</p>

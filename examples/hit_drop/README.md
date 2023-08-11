# Droplet in Homogeneous-Isotropic-Turbulence

- The density and viscosity ratios are fixed to 1. This allows us to use a FFT-based Poisson solver (adding a dependency on FFTW)
- The specific quantities to set in the `input` file are:
    - `Droplet diameter`, `position`, and `injection time` (which is expressed in terms of eddy-turnover times)
    - The `Weber number`, defined as in [Risso, F. & Fabre, J. 1998 Oscillations and breakup of a bubble immersed in a turbulent field. J. Fluid Mech. 372, 323–355](https://doi.org/10.1017/S0022112098002705)
$$\mathrm{We} = \dfrac{2 \rho \epsilon^{2/3} d^{5/3}}{\sigma}$$

    - The `Reynolds number` *(optional)* based on the Taylor microscale $\lambda$, defined as
$$\mathrm{Re_\lambda}=\dfrac{\mathrm{U_rms}\lambda}{\nu}$$

        - If not specified, $\mathrm{Re}_\lambda$ will be chosen as the maximum value allowed by the number of grid points.
    - The `Number of tracers` *(optional)*, specifying how many fluid tracers we want to be seeded on the droplet's surface.

<p align="center"><img src="./snap_t2.png"/></p>
<p align="center">Snaphost of the simulation at $t = 2\tau_L$ (after two eddy-turnover times).</p>

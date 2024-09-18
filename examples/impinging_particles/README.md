# Particle bed impinging into a pool

- This case simulates the impingment of a dense particle bed onto a water pool, following the experiments by [Saingier G, Sauret A, Jop P. Falling jet of dry granular material in water. Journal of Fluid Mechanics. 2021;916:A34.](https//doi.org/10.1017/jfm.2021.131)
- Particles are injected from a precomputed packed bed that is sorted prior to injection.
- The bed can be formed by running the `bedmaker` case using the `input.bedmaker` file as input.
- The case is dimensional.
- The case uses its own specialized `lpt_class.f90` file.


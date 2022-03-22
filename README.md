# NGA2

NGA2 is a high performance computing research library that provides a variety of finite volumes/finite difference solvers for typical fluid-related partial differential equations including:
- incompressible Navier-Stokes
- low-Mach number variable-density Navier-Stokes
- two-phase incompressible Navier-Stokes
- two-phase compressible Navier-Stokes
- constant and variable density scalar transport
- phasic volume fraction for volume-of-fluid methods
- overset and multi-block meshes via parallel coupler
- Lagrangian particle tracking
- large-eddy simulation models


Currently, NGA2 supports cartesian meshes only. However, because it is object-oriented, multiple meshes can be used simultaneously and interactively.

Future developments will focus on providing support for:
- cylindrical meshes
- dynamic remeshing
- basic unstructured meshes
- chemical kinetics and combustion models
- immersed boundaries

NGA2 is open-sourced under the [MIT license](./LICENSE).

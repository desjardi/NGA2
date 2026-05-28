# NGA2 Coding Rules

This document specifies coding rules and standards for the NGA2 codebase. These rules should be followed consistently across all modules.

## Reference Implementation

- **Use AMReX codebase as reference** (`~/Repositories/amrex`)
  - When implementing features similar to AMReX, check AMReX's implementation first
  - Compare nga2's implementation with AMReX's equivalent code
  - Do not invent solutions; verify against AMReX when available
  - AMReX is the authoritative reference for distributed data structures and communication patterns

## Documentation

- **DO NOT place documentation files in the repository**
  - Do not create `.md` files, README files, or other documentation in the repo unless explicitly requested
  - If you need to create documentation for reference, place it in `~/.cursor/docs/` directory
  - Only create documentation files in the repo when the user explicitly requests them

## MPI Commands

- **Use `mpiexec -n` instead of `mpirun -np`**
  - Example: `mpiexec -n 4 ./executable`
  - Not: `mpirun -np 4 ./executable`

## Real Number Precision

- **Use `real(WP)` instead of `real(8)`**
  - `WP` is provided by `use precision, only: WP`
  - Always include `use precision, only: WP` in modules that use real numbers

- **Use `_WP` suffix for all real literals**
  - Example: `1.0_WP`, `0.9_WP`, `0.0_WP`
  - Not: `1.0d0`, `0.9d0`, `0.0d0`

- **Use `real(..., WP)` for type conversions**
  - Example: `real(total_load, WP) / real(nprocs * max_load, WP)`
  - Not: `dble(total_load) / dble(nprocs * max_load)`

## Logical Operators

- **Prefer `.eq.` and `.ne.` over `==` and `!=`**
  - Example: `if (a .eq. b)` or `if (a .ne. b)`
  - Note: Both forms are supported via operator overloading, but `.eq.`/`.ne.` is preferred

## Indentation

- **Use 3 spaces for indentation**
  - Not tabs, not 2 spaces, not 4 spaces
  - Consistent 3-space indentation throughout

## Object-Oriented Design

- **Avoid constructors for derived types**
  - Use `initialize` type-bound procedures instead
  - Constructors can cause implicit copies and are "too dangerous"
  - Example: `call obj%initialize(...)` not `obj = type_name(...)`

- **All objects should have `initialize` and `finalize` type-bound procedures**
  - Users should not have to guess what needs to be linked or named
  - `initialize` should handle all setup (pointers, names, etc.)
  - `finalize` should handle all cleanup
  - Example: `call vtk%initialize(bg=bg, name='output')` not `vtk%bg => bg; vtk%name = 'output'`

## Parameter Configuration

- **Examples control parameter handling**
  - Modules provide default values for configuration parameters
  - Users read parameters from input files (using `inputfile_class` or `param`)
  - Users then call setter methods to configure modules
  - Modules should NOT automatically read from parameter systems
  - This follows nga2's design where examples control parameter handling

## Module Organization

- **Place `use` statements close to where they're used**
  - Use local `use` statements within subroutines/functions when possible
  - Exception: Keep module-level `use` statements if the imported symbols are used many times throughout the module
  - This improves code readability by making dependencies explicit at the point of use

## Comments and Documentation

- **Remove pedantic comments**
  - Avoid excessive mentions of AMReX where it's obvious
  - Keep comments focused and useful
  - Document design decisions, not obvious implementations

## Code Structure

- **Use block statements to encapsulate logical units in multistep routines**
  - Blocks should be used for **separate sequential steps** within a complex routine
  - Each block should contain variables declared **close to where they're used**, not all at the top
  - **DO NOT** wrap an entire routine in a single block - this defeats the purpose
  - Blocks help organize complex routines and improve readability by grouping related code and variables
  - Example:
    ```fortran
    subroutine complex_routine(this)
       implicit none
       class(mytype), intent(inout) :: this
       
       ! Step 1: Initialize
       init_block: block
          integer :: i, j
          ! ... initialization code using i, j ...
       end block init_block
       
       ! Step 2: Process data
       process_block: block
          real(WP) :: value
          ! ... processing code using value ...
       end block process_block
    end subroutine complex_routine
    ```
  - **Bad example** (DO NOT do this):
    ```fortran
    subroutine complex_routine(this)
       implicit none
       class(mytype), intent(inout) :: this
       
       everything_block: block
          integer :: i, j
          real(WP) :: value
          ! ... all code here ...
       end block everything_block
    end subroutine complex_routine
    ```

## Function/Subroutine Calls

- **Use explicit keyword arguments for input parameters**
  - Especially important for lengthy routine calls with many parameters
  - Improves readability and reduces errors
  - Example: `call bx1%initialize(lo=[0,0,0], hi=[31,31,31], itype=[CELL,CELL,CELL])`
  - Not: `call bx1%initialize([0,0,0], [31,31,31], [CELL,CELL,CELL])`

- **Data-first argument ordering for APIs**
  - When designing APIs, order arguments: required data first, then options, then optional name
  - Example: `call ens%add_scalar(data=pressure, comp=1, name="P")`
  - Example: `call amr%register(data=velocity, ncomp=3, ng=1, name="velocity")`
  - This makes the most important argument (the data) prominent

- **Use `target` attribute for objects stored by pointer**
  - When an object will be stored by pointer in a registry or collection, declare it with `target`
  - Example: `type(amrdata), target :: velocity`
  - The receiving procedure should use `class(...), target, intent(in)` for the argument

## Testing

- **Use the examples directory for tests**
  - Test cases should live in `examples/` subdirectories
  - Use the standard testing environment (not external test frameworks)
  - Tests should be comprehensive and cover all features
  - Run parallel tests with: `mpiexec -n 4 ./executable`

## AMR Implementation Notes (amrbase)

These are design choices in the current implementation, not strict rules:

- **Registry pattern**: Fields are registered with `amr%register(data=field, ...)` and automatically allocated/reallocated during regrid callbacks

- **Singleton pointer for callback routing**: AMReX callbacks don't pass user context, so `current_amrgrid` and `current_amrdata` module-level pointers are used to route callbacks to the correct object. This is a workaround for AMReX's C-style callback interface.

- **Boundary condition arrays**: `lo_bc(3,ncomp)` and `hi_bc(3,ncomp)` store per-direction, per-component BCs with defaults based on periodicity

## File Organization

- **Follow existing patterns**
  - Module files in appropriate `src/` subdirectories
  - Add new modules to `Make.package` files
  - Keep related functionality together

---

## Notes

- These rules are living guidelines and may be updated as the codebase evolves
- When in doubt, follow existing patterns in the codebase
- Consistency is more important than perfection - follow these rules even if alternatives exist
- Always check AMReX implementation when implementing similar features


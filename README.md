# OceanCutCells.jl

A Julia package for simulating ocean flows with irregular bathymetry using a cut-cell method and cell merging within the [Oceananigans.jl](https://clima.github.io/Oceananigans.jl/stable/) framework.

This package demonstrates a basic implementation of a cut-cell approach to represent complex seafloor topography on a rectilinear grid, including a simple cell merging strategy to improve numerical stability for small cut cells. Custom forcings are used to implement physics (advection, diffusion, pressure gradient, bottom drag) on this modified geometry.

## Installation

First, make sure you have [Julia installed](https://julialang.org/downloads/).

1.  Clone this repository:
    ```bash
    git clone https://github.com/your-username/OceanCutCells.jl.git # Replace with the actual repo URL
    cd OceanCutCells.jl
    ```
2.  Start Julia in the package directory and instantiate the environment:
    ```julia
    julia --project=.
    ```
    Once in the Julia REPL, enter package mode by typing `]` and run:
    ```julia
    pkg> instantiate
    ```
    This will download and install all dependencies listed in `Project.toml`.

## Running the Simulation

The `scripts/run_simulation.jl` file provides an example of how to set up and run a simulation using the `OceanCutCells.jl` package.

1.  Navigate to the package directory in your terminal.
2.  Run the script using Julia:
    ```bash
    julia --project=. scripts/run_simulation.jl
    ```

This script will:
*   Define simulation parameters.
*   Generate the grid and compute initial cut-cell geometry.
*   Perform the cell merging step.
*   Build the `Oceananigans.jl` model with custom forcings and geometry parameters.
*   Run the simulation for a specified duration.
*   Save the simulation output (fields and geometry) to a JLD2 file (`cut_cell_merged_simulation.jld2`).
*   (Optionally, if plotting is enabled in the script) Generate an animation (`CutCellTest_merged_pkg.gif`) from the saved output.

## Project Structure

The codebase is organized into a standard Julia package structure:

*   `src/`: Contains the core source code, split into modular files:
    *   `OceanCutCells.jl`: Main module file, includes other files and exports public functions/types.
    *   `Parameters.jl`: Defines the `CutCellParameters` struct holding all simulation and geometry information.
    *   `Geometry.jl`: Contains functions for bathymetry, computing initial geometry (`hFac` fields, wet areas/volumes), and the cell merging algorithm.
    *   `Forcings.jl`: Implements the custom physics forcings (advection, diffusion, PGF, bottom drag) that are "cut-cell aware".
    *   `Diagnostics.jl`: Implements custom diagnostics, specifically the continuity-based vertical velocity (`w`).
    *   `ModelSetup.jl`: Contains the `build_cut_cell_model` function that orchestrates the setup of the `Oceananigans.jl` model using the custom components.
*   `test/`: Contains unit tests.
    *   `runtests.jl`: The main test file that runs other test suites.
    *   `test_geometry.jl`: Tests specifically for the geometry computation and merging logic.
*   `scripts/`: Contains example scripts.
    *   `run_simulation.jl`: The main script to execute a simulation.
*   `Project.toml`: Specifies package dependencies.
*   `README.md`: This file.

## Customization and Development

*   **Parameters:** Adjust simulation parameters by modifying the `raw_params` struct in `scripts/run_simulation.jl`.
*   **Bathymetry:** Modify the `bathymetry_profile` and `bathymetry_derivative` functions in `src/Geometry.jl` to implement different seafloor shapes.
*   **Physics:** Extend or modify the custom forcing functions in `src/Forcings.jl` to change the active physics or their implementation details.
*   **Testing:** Run the test suite from the Julia REPL in the package directory using `] test`. Adding more tests (e.g., for forcings, model setup) is highly encouraged for development.

## Credits

This package is based on initial code and concepts for cut-cell implementation in Oceananigans.jl.

## License

[Consider adding a LICENSE.md file to specify the terms under which this code can be used and distributed, e.g., MIT License.]

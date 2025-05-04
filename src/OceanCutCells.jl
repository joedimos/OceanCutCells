module OceanCutCells

using Oceananigans
using Oceananigans.Units
using Oceananigans: Grid, VerticallyStretchedRectilinearGrid, architecture
using Oceananigans.Grids: xnodes, ynodes, znodes, topology, with_halo
using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel
using Oceananigans.TurbulenceClosures: nothing_closure
using Oceananigans.Advection: nothing_advection
using Oceananigans.Fields: ZeroField, Face, Center, ConstantField # Added ConstantField
using Oceananigans.Forcings: CustomForcing
using Oceananigans.BuoyancyModels: buoyancy_perturbation, SeawaterBuoyancy, LinearEquationOfState
using Oceananigans.Diagnostics: average, create_diagnostic_field, Diagnostic # Added Diagnostic
using Oceananigans.BoundaryConditions: fill_halo_regions, FluxBoundaryCondition, OpenBoundaryCondition # Added OpenBoundaryCondition
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, AbstractGridFittedBottom # Export IBG types and parent
using Oceananigans.Fields: Field # Export Field

using Printf # Moved from the script
using JLD2 # Moved from the script


include("Parameters.jl")
include("Geometry.jl") # Now includes CutCellBottom definition
include("Forcings.jl")
include("Diagnostics.jl")
include("ModelSetup.jl")


export
    # Parameters
    CutCellParameters, # Now contains only physical parameters
    # Geometry
    bathymetry_profile,
    bathymetry_derivative, # Keep in case needed
    CutCellBottom, # Export the new immersed boundary type
   
    
    add_cut_cell_pressure_gradient_force!,
    add_cut_cell_advection!, # Generic advection
    add_cut_cell_diffusion!, # Generic diffusion
    add_cut_cell_vertical_advection_u!,
    add_cut_cell_vertical_diffusion_u!,
    add_cut_cell_bottom_drag!,

    # Diagnostics
    diagnose_cut_cell_w!,

    # Model Setup
    build_cut_cell_model,

    # Re-export useful Oceananigans types/functions/macros
    VerticallyStretchedRectilinearGrid,
    ImmersedBoundaryGrid, # Export ImmersedBoundaryGrid
    Simulation,
    JLD2OutputWriter,
    IterationInterval,
    Callback,
    prettytime,
    time,
    iteration,
    set!,
    maximum,
    abs,
    FieldBoundaryConditions,
    FluxBoundaryCondition,
    OpenBoundaryCondition,
    Field, # Export Field for creating bathymetry etc.
    Center, # Export grid locations
    Face,
    @info # Export @info for logging

end # module OceanCutCells

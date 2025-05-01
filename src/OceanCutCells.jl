    module OceanCutCells

    using Oceananigans
    using Oceananigans.Units
    using Oceananigans: Grid, VerticallyStretchedRectilinearGrid, architecture
    using Oceananigans.Grids: xnodes, ynodes, znodes
    using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel
    using Oceananigans.TurbulenceClosures: HorizontalLeapfrogDiffusivity, ScalarDiffusivity, nothing_closure
    using Oceananigans.Advection: CenteredSecondOrder, UpwindBiasedThirdOrder, nothing_advection
    using Oceananigans.Fields: ZeroField, Face, Center, ConstantField # Added ConstantField
    using Oceananigans.Forcings: CustomForcing
    using Oceananigans.BuoyancyModels: buoyancy_perturbation, SeawaterBuoyancy, LinearEquationOfState
    using Oceananigans.Diagnostics: average, create_diagnostic_field, Diagnostic # Added Diagnostic
    using Oceananigans.BoundaryConditions: fill_halo_regions, FluxBoundaryCondition, OpenBoundaryCondition # Added OpenBoundaryCondition

    using Printf # Moved from the script
    using JLD2 # Moved from the script

    
    include("Parameters.jl")
    include("Geometry.jl")
    include("Forcings.jl")
    include("Diagnostics.jl")
    include("ModelSetup.jl")

    
    export
        # Parameters
        CutCellParameters,
        # Geometry
        bathymetry_profile,
        bathymetry_derivative,
        compute_initial_geometry,
        merge_cut_cells,
        # Forcings (users might need to define custom ones)
        add_cut_cell_pressure_gradient_force!,
        add_cut_cell_advection_T!,
        add_cut_cell_diffusion_T!,
        add_cut_cell_vertical_advection_u!,
        add_cut_cell_vertical_diffusion_u!,
        add_cut_cell_bottom_drag!,
        # Diagnostics
        diagnose_cut_cell_w!,
        # Model Setup
        build_cut_cell_model,
        # Re-export useful Oceananigans types/functions (optional, but convenient)
        VerticallyStretchedRectilinearGrid,
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
        OpenBoundaryCondition # Export OpenBC if used

    end # module OceanCutCells

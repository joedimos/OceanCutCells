    module ModelSetup

    using Oceananigans
    using Oceananigans.Units
    using Oceananigans.Grids: VerticallyStretchedRectilinearGrid, xnodes, znodes
    using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel
    using Oceananigans.TurbulenceClosures: nothing_closure
    using Oceananigans.Advection: nothing_advection
    using Oceananigans.Fields: Center, Face, ZeroField, ConstantField # Added ConstantField
    using Oceananigans.BuoyancyModels: SeawaterBuoyancy, LinearEquationOfState
    using Oceananigans.BoundaryConditions: FieldBoundaryConditions, FluxBoundaryCondition, OpenBoundaryCondition
    using Oceananigans.Diagnostics: create_diagnostic_field, Diagnostic

    using ..Parameters: CutCellParameters # Access the struct
    using ..Geometry: compute_initial_geometry, merge_cut_cells, bathymetry_profile, bathymetry_derivative # Need geometry functions
    using ..Forcings: u_forcings, T_forcings, S_forcings # Need forcing definitions
    using ..Diagnostics: diagnose_cut_cell_w! # Need diagnostic function

    export build_cut_cell_model # Export the main function

    function build_cut_cell_model(params::CutCellParameters)

        # --- Grid Generation ---
        # Non-uniform z-spacing: More resolution near the surface
        # This function definition should maybe live here or be passed in.
        # For now, let's define it here as it's specific to this grid setup.
        z_stretching = 1.3 # Factor for vertical stretching
        z_faces_func(k) = -params.total_depth * (k/params.nz)^z_stretching

        # Create a VerticallyStretchedRectilinearGrid
        grid = VerticallyStretchedRectilinearGrid(
            topology = (Bounded, Flat, Bounded),
            size = (params.nx, params.nz),
            x = params.x_domain,
            z_faces = z_faces_func.(0:params.nz),
            halo = (2, 2) # 2 in x, 2 in z
        )
        @info "Grid created: $(grid)"

        # --- Cut Cell Geometry: hFac fields and Wet Areas/Volumes ---
        @info "Computing initial cut cell geometry..."
        hFacC_initial, hFacW_initial, hFacS_initial,
        wet_cell_volume_initial, wet_face_area_x_initial, wet_face_area_z_initial =
            compute_initial_geometry(grid, bathymetry_profile) # Use the function from Geometry.jl

        # --- Implement Cell Merging ---
        # Use the function from Geometry.jl
        hFacC_merged, wet_cell_volume_merged, recip_wet_cell_volume_merged,
        hFacW_geom, hFacS_geom, wet_face_area_x_geom, wet_face_area_z_geom =
             merge_cut_cells(grid, hFacC_initial, wet_cell_volume_initial, bathymetry_derivative,
                             hFacW_initial, hFacS_initial, wet_face_area_x_initial, wet_face_area_z_initial)

        # Bundle Cut-Cell Parameters (using merged geometry where appropriate)
        # Pass the geometry fields back to the parameters struct.
        # Note: We are using the initial (non-merged) hFacW/S and face areas
        # as the merging logic only applied to C cells and volumes.
        cc_params_with_geom = CutCellParameters(
            hFacC = hFacC_merged,
            hFacW = hFacW_geom, # These come back from merge_cut_cells
            hFacS = hFacS_geom, # These come back from merge_cut_cells
            wet_cell_volume = wet_cell_volume_merged,
            recip_wet_cell_volume = recip_wet_cell_volume_merged,
            wet_face_area_x = wet_face_area_x_geom, # These come back from merge_cut_cells
            wet_face_area_z = wet_face_area_z_geom, # These come back from merge_cut_cells
            rho0 = params.rho0,
            g = params.g,
            Kh = params.Kh,
            Kv = params.Kv,
            Ah = params.Ah,
            Av = params.Av,
            bottom_drag_coeff = params.bottom_drag_coeff,
            alpha = params.alpha, # Added parameters
            beta = params.beta,
            T0 = params.T0,
            S0 = params.S0,
            nx = params.nx, # Added parameters
            nz = params.nz,
            total_width = params.total_width,
            total_depth = params.total_depth,
            x_domain = params.x_domain,
            z_domain = params.z_domain,
            dt = params.dt # Added parameters
        )

        # --- Model Setup ---
        # Set standard advection and closure to 'nothing' to use custom forcings exclusively
        advection = nothing_advection()
        closure = nothing_closure()

        buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=params.alpha, β=params.beta, T₀=params.T0, S₀=params.S0))

        # Boundary conditions: Solid wall (zero velocity, zero flux) on sides (x) and bottom (z) by default for Bounded.
        # Top (z=0) BC for u: Wind stress applied as a kinematic flux.
        wind_stress_tau_x = 0.05 # Should be a parameter
        wind_stress_flux = wind_stress_tau_x / params.rho0

        u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(wind_stress_flux))

        # No explicit BCs needed for T and S with custom diffusion forcings if flux = -K*grad is used and grad is handled at boundary.
        # For a passive tracer/heat, no flux at solid boundaries is the default zero gradient (∂T/∂n = 0).
        # Our custom diffusion handles this by calculating gradient only if both adjacent cells are wet.
        # For the top surface, need a BC for T/S. Let's assume zero flux (insulation).
        # The custom diffusion function `add_cut_cell_diffusion!` implicitly handles this for the top face (k=1, face k=1)
        # by requiring `hFacS[i,j,k-1]` (which is `hFacS[i,j,0]`, invalid index) and wet cells `k-1`, `k`.
        # Since `k-1` is not in the grid for k=1, the gradient calculation `is_wet_ikm1 && is_wet_ik` fails for the top face.
        # This correctly results in zero flux across the top face via the custom forcing logic.

        model = HydrostaticFreeSurfaceModel(
            grid = grid,
            advection = advection, # Standard advection is off
            coriolis = nothing, # Add Coriolis if needed
            buoyancy = buoyancy,
            closure = closure,     # Standard diffusion/viscosity is off
            boundary_conditions = (u=u_bcs,), # Apply wind stress BC to the u field
            tracers = (:T, :S),
            parameters = (; cut_cell_params = cc_params_with_geom), # Pass the parameters including geometry
            forcings = (u=u_forcings, T=T_forcings, S=S_forcings), # Apply all custom forcings
        )

        # --- Initial Conditions ---

        # Linear stratification
        T_init(x, y, z) = params.T0 - 10.0 * (z - grid.zᵃᵃᶠ[1]) / (grid.zᵃᵃᶠ[grid.Nz+1] - grid.zᵃᵃᶠ[1]) # Colder deeper
        S_init(x, y, z) = params.S0 + 0.5 * (z - grid.zᵃᵃᶠ[1]) / (grid.zᵃᵃᶠ[grid.Nz+1] - grid.zᵃᵃᶠ[1]) # Saltier deeper

        # Set initial temperature and salinity
        set!(model, T=T_init, S=S_init)

        # Initialize velocities to zero
        set!(model, u=0.0, v=0.0, w=0.0) # Model w field is updated by diagnostic

        # Apply masks to initial conditions based on the *merged* hFacC
        # Ensure tracers are zero in dry cells (including those merged away).
        @info "Masking initial T and S fields with merged hFacC..."
        @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
            if cc_params_with_geom.hFacC[i, j, k] < 1e-10 # Use the merged hFacC
                # Use direct data access with halo indices for speed if modifying many points
                model.tracers.T.data[i+grid.Hx, j+grid.Hy, k+grid.Hz] = 0.0
                model.tracers.S.data[i+grid.Hx, j+grid.Hy, k+grid.Hz] = 0.0
            end
        end
        # Velocities should also be zero in dry/solid regions.
        # u is at Face x. If hFacW=0, u should be 0. Mask using the non-merged hFacW.
        @info "Masking initial u field with hFacW..."
        @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
             if cc_params_with_geom.hFacW[i, j, k] < 1e-10
                 model.velocities.u.data[i+grid.Hx-1, j+grid.Hy, k+grid.Hz] = 0.0
             end
        end
        # w is at Face z. If hFacS=0, w should be 0. Mask using the non-merged hFacS.
        @info "Masking initial w field with hFacS..."
        @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz+1
             if cc_params_with_geom.hFacS[i, j, k] < 1e-10
                 # The model w field is not computed directly, it's a diagnostic output.
                 # We don't need to mask the *model's* internal w (which doesn't exist or is always zero),
                 # but we need to mask the *diagnostic* w field *after* it's computed.
                 # The diagnostic function itself handles setting w=0 where hFacS=0.
                 # So no need to mask the initial model.velocities.w (which is likely a ZeroField or ConstantField(0)).
             end
        end

        # --- Custom Diagnostics ---
        # The w diagnostic calculates w *after* the state update using u, then stores it in cut_cell_w.
        cut_cell_w = create_diagnostic_field(model.velocities.w) # Create a field to store the diagnosed w
        w_diagnostic = Diagnostic(diagnose_cut_cell_w!, model, field_dependencies = (model.velocities.u,), field_outputs = (w=cut_cell_w,))

        # Add the w diagnostic to the model
        model.diagnostics[:cut_cell_w] = w_diagnostic # Store the diagnostic object itself

        return model, cut_cell_w # Return the model and the diagnostic field
    end

    end # module ModelSetup

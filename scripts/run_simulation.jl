    using Oceananigans
    using Oceananigans.Units
    using JLD2
    using CairoMakie # Needed for plotting later
    using Printf

    
    using OceanCutCells
    using OceanCutCells.Parameters
    using OceanCutCells.ModelSetup
    using OceanCutCells.Geometry # Needed for plotting grid coords and bathymetry

    # --- Simulation Parameters ---
    # Define simulation-specific parameters here, maybe overriding defaults
    params = CutCellParameters( # Fill in all fields from the struct
        nx = 40,
        nz = 20,
        total_width = 40000.0,
        total_depth = 3000.0,
        x_domain = (0.0, 40000.0), # Ensure Float64
        z_domain = (-3000.0, 0.0), # Ensure Float64

        dt = 100.0,           # Time step
        Kh = 5.0,
        Kv = 1e-4,
        Ah = 20.0,
        Av = 1e-3,
        rho0 = 1025.0,
        g = 9.81,
        alpha = 2e-4,
        beta = 7.4e-4,
        T0 = 15.0,
        S0 = 35.5,
        bottom_drag_coeff = 2e-3,
        # These fields are filled during geometry setup, just dummy values here
        hFacC = nothing, hFacW = nothing, hFacS = nothing,
        wet_cell_volume = nothing, recip_wet_cell_volume = nothing,
        wet_face_area_x = nothing, wet_face_area_z = nothing,
    )

    stop_time = 1000 * params.dt # Run for some time
    output_interval = 50 * params.dt # Save output every X time units
    filename = "cut_cell_merged_output_pkg.jld2"


    # --- Build Model ---
    # build_cut_cell_model will handle grid creation, geometry, forcings, diagnostics, ICs
    @info "Building model..."
    model, diagnosed_w_field = build_cut_cell_model(params)
    @info "Model built."

    # --- Simulation ---
    simulation = Simulation(model, Δt=params.dt, stop_time=stop_time)

    # --- Output Writer ---
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model,
        (T=model.tracers.T, u=model.velocities.u, S=model.tracers.S, p=model.pressure, # p is standard output
         hFacC=model.parameters.cut_cell_params.hFacC, # Access parameters via model
         wet_cell_volume=model.parameters.cut_cell_params.wet_cell_volume,
         wet_face_area_x=model.parameters.cut_cell_params.wet_face_area_x,
         wet_face_area_z=model.parameters.cut_cell_params.wet_face_area_z,
        ),
        other_fields = (; w=diagnosed_w_field), # Output our custom 'w' field
        filename = filename,
        schedule = TimeInterval(output_interval), # Use TimeInterval
        overwrite_existing = true
    )

    # --- Progress Reporting ---
    progress(sim) = @printf("i: %d, t: %s, wall time: %s, max|u|: %.2e, max|w|: %.2e, max|T|: %.2f\n",
                             iteration(sim), prettytime(time(sim)), prettytime(sim.run_wall_time),
                             maximum(abs, sim.model.velocities.u), maximum(abs, sim.diagnostics[:cut_cell_w].field), maximum(abs, sim.model.tracers.T)) # Access w via diagnostic object

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(50)) # Report every 50 iterations

    # --- Run the simulation ---
    println("Starting simulation with cell merging...")
    run!(simulation)
    println("Simulation complete.")

    # --- Analysis and Plotting (Optional - could go in a separate script) ---

    println("Generating animation...")

    file = jldopen(filename)
    iterations = parse.(Int, keys(file["timeseries/t"])) # Use keys for iterations
    grid = file["grid"] # Load grid from file

    # Extract coordinate arrays
    xC = xnodes(Center(), grid)[:] # Simplified access if grid is loaded
    zC = znodes(Center(), grid)[:]
    xF = xnodes(Face(), grid)[:]
    # Need z locations corresponding to u and w
    zCu = znodes(Face(), grid, Center())[:] # z at u (C)
    zCw = znodes(Center(), grid, Face())[:] # z at w (F)

    # Get bathymetry for plotting
    x_bath = xnodes(Center(), grid)[:,1,1] # x locations at cell centers
    bathymetry_depth = bathymetry_profile.(x_bath) # Use the function from Geometry module

    # --- Create animation ---
    # Use the loaded grid for fields
    anim = @animate for i in iterations
        @info "Plotting iteration $i..."
        T = file["timeseries/T/$i"][:, 1, :]
        u = file["timeseries/u/$i"][:, 1, :]
        w = file["timeseries/w/$i"][:, 1, :] # Get the custom diagnosed w
        hFacC = file["timeseries/hFacC/$i"][:, 1, :]
        p = file["timeseries/p/$i"][:, 1, :] # Use 'p' saved from model.pressure
        wet_vol = file["timeseries/wet_cell_volume/$i"][:, 1, :]

        # Plotting
        fig = Figure(resolution = (1200, 2100)) # Taller figure for 6 plots
        ax_T = Axis(fig[1, 1], title = "Temperature (°C)")
        ax_u = Axis(fig[1, 2], title = "Velocity u (m/s)")
        ax_w = Axis(fig[2, 1], title = "Velocity w (m/s)")
        ax_p = Axis(fig[2, 2], title = "Pressure Anomaly (Pa)")
        ax_hFacC = Axis(fig[3, 1], title = "Merged hFacC") # Show merged hFacC
        ax_wetvol = Axis(fig[3, 2], title = "Merged Wet Volume")

        # Use `heatmap!` which handles the grid and data indexing correctly
        hm_T = heatmap!(ax_T, xC, zC, T, colormap = :RdBu_r)
        Colorbar(fig[1, 1][2, 1], hm_T)

        hm_u = heatmap!(ax_u, xF, zCu, u, colormap = :RdBu_r, colorrange=(-0.1, 0.1))
        Colorbar(fig[1, 2][2, 1], hm_u)

        hm_w = heatmap!(ax_w, xC, zCw, w, colormap = :RdBu_r, colorrange=(-5e-4, 5e-4))
        Colorbar(fig[2, 1][2, 1], hm_w)

        hm_p = heatmap!(ax_p, xC, zC, p, colormap = :viridis) # Use xC, zC for pressure (at C)
        Colorbar(fig[2, 2][2, 1], hm_p)

        hm_hFacC = heatmap!(ax_hFacC, xC, zC, hFacC, colormap = :Blues, colorrange=(0, 1))
        Colorbar(fig[3, 1][2, 1], hm_hFacC)

        hm_wetvol = heatmap!(ax_wetvol, xC, zC, wet_vol, colormap = :Plasma)
        Colorbar(fig[3, 2][2, 1], hm_wetvol)

        # Overlay bathymetry profile on the T and hFacC plots
        lines!(ax_T, x_bath, -bathymetry_depth, color=:black, linewidth=2)
        lines!(ax_hFacC, x_bath, -bathymetry_depth, color=:black, linewidth=2)
        lines!(ax_wetvol, x_bath, -bathymetry_depth, color=:black, linewidth=2)


        # Add a title above the plots indicating the time
        Label(fig[0, :], text = "Cut Cell Test with Merging (Flux-Based Forcings): " * @sprintf("%.2f", file["timeseries/t/$i"][end] / day) * " days",
              fontsize = 30)

        # Set axis labels and limits
        for ax in [ax_T, ax_u, ax_w, ax_p, ax_hFacC, ax_wetvol]
            ax.xlabel = "x (m)"
            ax.ylabel = "z (m)"
            ylims!(ax, -params.total_depth, 0)
            xlims!(ax, 0, params.total_width)
        end

        fig
    end

    save("CutCellTest_merged_pkg.gif", anim, fps = 10)
    close(file)
    println("Animation saved in CutCellTest_merged_pkg.gif")

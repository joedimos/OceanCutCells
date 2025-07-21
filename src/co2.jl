using Oceananigans
using Oceananigans.Units: minutes, hour
using Oceananigans.Advection: UpwindBiasedFifthOrder
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Statistics
using Printf
using JLD2
using Plots

# ----------------------
# 1. Grid and Immersed Boundary (DAC Collectors)
# ----------------------

# A regular rectilinear grid that underlies the immersed boundary grid
nx, ny = 100, 50
Lx, Ly = 100.0, 50.0
underlying_grid = RectilinearGrid(size=(nx, ny, 1),
                                  extent=(Lx, Ly, 1.0),
                                  topology=(Bounded, Bounded, Flat))

# Define the DAC collectors as circular immersed boundaries
const DAC_LOCATIONS = [(50.0, 25.0), (70.0, 10.0)]
const DAC_RADIUS = 2.0

# This function returns `true` if a grid point is inside a DAC collector, marking it as "land"
@inline function is_in_dac(x, y, z)
    for (xc, yc) in DAC_LOCATIONS
        if (x - xc)^2 + (y - yc)^2 < DAC_RADIUS^2
            return true
        end
    end
    return false
end

# The ImmersedBoundaryGrid wraps the underlying grid with the DAC geometry
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(is_in_dac))

# ----------------------
# 2. Model components
# ----------------------
# Constant rightward wind (velocity field)
u_vel(x, y, z, t) = 1.0
v_vel(x, y, z, t) = 0.0
w_vel(x, y, z, t) = 0.0
velocities = (u=u_vel, v=v_vel, w=w_vel)

# DAC sink parameters
const CO2_TARGET = 400.0        # ppm target threshold
const GAIN = 0.2                # Proportional gain of removal
const MAX_REMOVAL = 2.0 / 60    # ppm per second max removal rate (converted from per minute)
const ENERGY_PER_PPM = 5.0      # Arbitrary energy units per ppm removed

# This forcing function implements the feedback-controlled CO₂ sink.
# It is only active at the boundary of the immersed DAC collectors.
@inline function dac_co2_removal(i, j, k, grid, clock, fields)
    # Check if the current cell (i, j, k) is an immersed boundary cell
    # The condition `grid.immersed_boundary.bottom_height[i, j] > 0` is true for cells over an immersed boundary in a flat model.
    if grid.immersed_boundary.bottom_height[i, j] > 0
        c = @inbounds fields.c[i, j, k]
        if c > CO2_TARGET
            # Calculate sink based on proportional control
            sink = -GAIN * (c - CO2_TARGET)
            # Return the sink, capped by the maximum removal rate
            return max(sink, -MAX_REMOVAL)
        end
    end
    return 0.0
end

# The forcing is correctly applied to the tracer 'c' in `sources`
c_forcing = Forcing(dac_co2_removal, discrete_form=true)

# ----------------------
# 3. Initialize Oceananigans model
# ----------------------
model = NonhydrostaticModel(grid=grid,
                            advection=UpwindBiasedFifthOrder(),
                            timestepper=:RungeKutta3,
                            tracers=:c,
                            velocities=velocities,
                            sources=(c=c_forcing,))

# ----------------------
# 4. Set initial CO₂ field
# ----------------------
# Elevated CO₂ on the left side, representing an incoming plume
set!(model, c = (x, y, z) -> x < 20 ? 420.0 : 400.0)

# ----------------------
# 5. Simulation setup
# ----------------------
simulation = Simulation(model, Δt=1.0, stop_time=100minutes)

# Create a field to store the calculated sink term for energy tracking
S_field = Field{Center, Center, Center}(model.grid)
compute!(S_field, model.sources.c.func)
# This callback calculates and stores the total energy consumed by the DACs
global_energy = Float64[]
function energy_callback(sim)
    compute!(S_field, model.sources.c.func)
    # The total removal is the integral of the sink term over the volume
    # Since sink is in ppm/s, we get total ppm/s removed.
    # Volume of a cell is Lx/nx * Ly/ny * 1.0
    cell_volume = (Lx/nx) * (Ly/ny) * 1.0
    total_removal_rate = abs(sum(S_field) * cell_volume)
    # Energy is the removal rate multiplied by the energy cost
    energy = ENERGY_PER_PPM * total_removal_rate
    push!(global_energy, energy)
end

# Output writer for saving CO₂ data
simulation.output_writers[:co2_output] = JLD2OutputWriter(model, model.tracers,
                                                          schedule=TimeInterval(10minutes),
                                                          prefix="co2_cutcell_feedback",
                                                          force=true)
# Callback for printing progress
simulation.callbacks[:progress] = Callback(sim ->
    @info @sprintf("Time = %s | Mean CO₂ = %.2f ppm",
                   prettytime(sim.clock.time),
                   mean(sim.model.tracers.c)),
                   TimeInterval(10minutes))

# Callback for logging energy usage at each time step
simulation.callbacks[:energy_log] = Callback(energy_callback, IterationInterval(1))

# ----------------------
# 6. Run the simulation
# ----------------------
run!(simulation)
@info "Simulation complete."

# ----------------------
# 7. Visualization
# ----------------------
# Load the results from the output file
filepath = "co2_cutcell_feedback.jld2"
file = jldopen(filepath)
iterations = keys(file["timeseries/c"])
final_iter = iterations[end]
c_final = file["timeseries/c/$final_iter"]
close(file)

# Get grid coordinates for plotting
xc, yc, zc = nodes(model.tracers.c)

# Plot the final CO₂ concentration
heatmap(xc, yc, c_final[:, :, 1]',
        aspect_ratio=:equal,
        xlabel="x (m)",
        ylabel="y (m)",
        title="CO₂ Concentration at t = $(prettytime(simulation.stop_time)) (ppm)",
        c=:thermal,
        clims=(398, 421))

# Plot the energy consumption over time
plot(global_energy,
     xlabel="Timestep",
     ylabel="Energy Rate (units/s)",
     title="DAC Energy Consumption Rate",
     legend=false)
using Oceananigans
using Oceananigans.Units: minutes, hour
using Oceananigans.Advection: UpwindBiasedFifthOrder
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary

using Statistics
using Printf
using JLD2
using Plots

----------------------------------------------------------------------
#         DAC Collectors
----------------------------------------------------------------------

# A regular rectilinear grid provides the foundation for our simulation domain.
# This underlying grid does not yet know about the DAC collectors.
nx, ny = 200, 100
Lx, Ly = 200.0, 100.0 # Domain size in meters
underlying_grid = RectilinearGrid(size=(nx, ny, 1),
                                  extent=(Lx, Ly, 1.0),
                                  topology=(Bounded, Bounded, Flat))

# --- Define the geometry of the DAC collectors ---
# Collectors at specific (x, y) locations with a defined radius.
const DAC_LOCATIONS = [(60.0, 50.0), (100.0, 25.0), (100.0, 75.0)]
const DAC_RADIUS = 5.0

# This function defines the shape of the immersed boundary. It returns `true`
# if a given coordinate (x, y, z) is inside any of the DAC collectors,
# effectively marking that region as a solid "land" area that fluid cannot pass through.
@inline function is_in_dac(x, y, z)
    for (xc, yc) in DAC_LOCATIONS
        if (x - xc)^2 + (y - yc)^2 < DAC_RADIUS^2
            return true
        end
    end
    return false
end

# The ImmersedBoundaryGrid combines the `underlying_grid` with the DAC geometry.

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(is_in_dac))

# ----------------------------------------------------------------------
#          
       PHYSICS AND FORCING----------------------------------------------------------------------

# --- Define a constant, steady wind from left to right ---

# This represents a constant background wind flowing over the DAC facility.
u_vel(x, y, z, t) = 1.0    # 1.0 m/s in the x-direction
v_vel(x, y, z, t) = 0.0
w_vel(x, y, z, t) = 0.0
velocities = (u=u_vel, v=v_vel, w=w_vel)

# --- Define the CO₂ removal logic for the DAC collectors ---
const CO2_BACKGROUND = 400.0       # Background CO₂ concentration (ppm)
const CO2_PLUME = 420.0            # Plume CO₂ concentration (ppm)
const CO2_TARGET = 400.0           # Target threshold for DAC activation (ppm)
const GAIN = 0.1                   # Proportional gain for the feedback control
const MAX_REMOVAL = 2.0 / 60       # Max removal rate (ppm/s), converted from per minute
const ENERGY_PER_PPM = 5.0         # Arbitrary energy cost per ppm of CO₂ removed

# This forcing function implements the feedback-controlled CO₂ sink.
# It is only active at the grid cells that form the boundary of the DAC collectors.
@inline function dac_co2_removal(i, j, k, grid, clock, fields)
    # This condition identifies a cell adjacent to an immersed boundary.
    # In our flat (z=1) grid, `bottom_height` is non-zero for cells "on" the boundary.
    if grid.immersed_boundary.bottom_height[i, j] > 0
        c = @inbounds fields.c[i, j, k] # Get current CO₂ concentration

        # Only activate the sink if the concentration exceeds the target
        if c > CO2_TARGET
            # The removal rate is proportional to the excess CO₂.
            # This is a simple "proportional control" system.
            sink = -GAIN * (c - CO2_TARGET)

            # Cap the removal rate to the maximum physical capacity of the DAC.
            # `max` is used because `sink` is negative.
            return max(sink, -MAX_REMOVAL)
        end
    end

    # If not on a boundary or if CO₂ is below target, do nothing.
    return 0.0
end

# Create the forcing object to be used in the model
c_forcing = Forcing(dac_co2_removal, discrete_form=true)

# ----------------------------------------------------------------------
#                  SECTION 3: MODEL INSTANTIATION AND SETUP
# ----------------------------------------------------------------------

# Advection: UpwindBiasedFifthOrder is a good choice for minimizing numerical diffusion
#            while accurately transporting the CO₂ plume.
# Timestepper: RungeKutta3 is a standard, stable, and accurate time-stepping scheme.
model = NonhydrostaticModel(grid=grid,
                            advection=UpwindBiasedFifthOrder(),
                            timestepper=:RungeKutta3,
                            tracers=:c, # 'c' represents CO₂ concentration
                            velocities=velocities,
                            forcing=(c=c_forcing,))

# --- Set the initial CO₂ field ---
# Plume of high-concentration CO₂
# is located on the left side of the domain, ready to be blown into the DACs.
set!(model, c = (x, y, z) -> x < 30 ? CO2_PLUME : CO2_BACKGROUND)

# ----------------------------------------------------------------------
#                     SECTION 4: SIMULATION SETUP
# ----------------------------------------------------------------------

simulation = Simulation(model, Δt=0.5, stop_time=200minutes)

# --- Set up callbacks for data collection and progress logging ---

# Create a field to store the calculated sink term for energy tracking
S_field = Field{Center, Center, Center}(model.grid)

# Arrays to store the time-series data for the energy plot
global_energy = Float64[]
simulation_times = Float64[]

# This callback calculates and stores the total energy consumed by the DACs at each step.
function energy_callback(sim)
    # forcing term `model.forcing.c'
    #  `S_field`. The original code `model.sources.c.func`
    
    compute!(S_field, model.forcing.c)

    # The total removal rate is the integral of the sink term over the whole volume.
    # The volume of a cell is constant in a RectilinearGrid.
    cell_volume = (Lx/nx) * (Ly/ny) * 1.0
    total_removal_rate = abs(sum(S_field) * cell_volume)

    # Energy consumption is the removal rate multiplied by the energy cost.
    energy = ENERGY_PER_PPM * total_removal_rate
    push!(global_energy, energy)
    push!(simulation_times, sim.clock.time)
end

# Output writer for saving the full 3D CO₂ field at regular intervals
simulation.output_writers[:co2_output] = JLD2OutputWriter(model, model.tracers,
                                                          schedule=TimeInterval(1minute),
                                                          prefix="co2_dac_feedback",
                                                          force=true)

# Callback for printing a progress message in the terminal
simulation.callbacks[:progress] = Callback(sim ->
    @info @sprintf("Time = %s | Mean CO₂ = %.2f ppm | Max CO₂ = %.2f ppm",
                   prettytime(sim.clock.time),
                   mean(sim.model.tracers.c),
                   maximum(sim.model.tracers.c)),
                   TimeInterval(10minutes))

# Add the energy calculation callback to the simulation, running at every iteration
simulation.callbacks[:energy_log] = Callback(energy_callback, IterationInterval(1))

# ----------------------------------------------------------------------
#                        SECTION 5: RUN SIMULATION
# ----------------------------------------------------------------------
@info "Simulation starting..."
run!(simulation)
@info "Simulation complete."

# ----------------------------------------------------------------------
#                SECTION 6: POST-PROCESSING & VISUALIZATION
# ----------------------------------------------------------------------
@info "Generating visualizations..."

# --- Load data from the saved output file ---
filepath = "co2_dac_feedback.jld2"
c_timeseries = FieldTimeSeries(filepath, "c")
times = c_timeseries.times

# Get grid coordinates for plotting
xc, yc, zc = nodes(c_timeseries[1])

# --- Plot 1: Final CO₂ Concentration with DAC locations ---
final_time_index = length(times)
c_final = c_timeseries[final_time_index]

# Function to draw a circle for plotting
function circle_shape(h, k, r)
    θ = range(0, 2π, length=100)
    return h .+ r*sin.(θ), k .+ r*cos.(θ)
end

heatmap(xc, yc, interior(c_final, :, :, 1)',
        aspect_ratio=:equal,
        xlabel="x (meters)",
        ylabel="y (meters)",
        title="Final CO₂ Concentration at t = $(prettytime(times[end]))",
        c=:thermal,
        clims=(CO2_BACKGROUND-2, CO2_PLUME+1))

# Overlay the DAC collector locations on the heatmap
for (xc_dac, yc_dac) in DAC_LOCATIONS
    plot!(circle_shape(xc_dac, yc_dac, DAC_RADIUS),
          seriestype=[:shape,], lw=2, c=:white, linecolor=:black,
          legend=false, fillalpha=0.2)
end
savefig("final_co2_concentration.png")


# --- Plot 2: Energy Consumption Rate Over Time ---
plot(simulation_times, global_energy,
     xlabel="Time (seconds)",
     ylabel="Energy Rate (units/s)",
     title="DAC Energy Consumption Rate Over Time",
     legend=false,
     linewidth=2)
savefig("dac_energy_consumption.png")


# --- Animation: CO₂ Field Evolution ---
anim = @animate for i in 1:length(times)
    @info "Animating frame $i of $(length(times))..."
    c_t = c_timeseries[i]
    heatmap(xc, yc, interior(c_t, :, :, 1)',
            aspect_ratio=:equal,
            xlabel="x (meters)",
            ylabel="y (meters)",
            title="CO₂ Concentration at t = $(prettytime(times[i])) (ppm)",
            c=:thermal,
            clims=(CO2_BACKGROUND-2, CO2_PLUME+1))

    # Overlay the DAC collector locations
    for (xc_dac, yc_dac) in DAC_LOCATIONS
        plot!(circle_shape(xc_dac, yc_dac, DAC_RADIUS),
              seriestype=[:shape,], lw=2, c=:white, linecolor=:black,
              legend=false, fillalpha=0.2)
    end
end

gif(anim, "co2_dac_evolution.gif", fps = 15)

@info "Visualization complete. Files saved: final_co2_concentration.png, dac_energy_consumption.png, co2_dac_evolution.gif"
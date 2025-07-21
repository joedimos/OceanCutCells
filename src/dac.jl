using Oceananigans
using Oceananigans.Units: minutes, hour
using Oceananigans.Advection: UpwindBiasedFifthOrder
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.OutputWriters
using Statistics
using Printf
using JLD2
using Plots

# ----------------------------------------------------------------------------------
# DOMAIN AND GRID
# ----------------------------------------------------------------------------------

nx, ny = 200, 100
Lx, Ly, Lz = 200.0, 100.0, 1.0
underlying_grid = RectilinearGrid(size=(nx, ny, 1),
extent=(Lx, Ly, Lz),
topology=(Bounded, Bounded, Flat))

# DAC Collectors: locations and radius
const DAC_LOCATIONS = [(60.0, 50.0), (100.0, 25.0), (100.0, 75.0)]
const DAC_RADIUS = 5.0

@inline function is_in_dac(x, y, z)
for (xc, yc) in DAC_LOCATIONS
if (x - xc)^2 + (y - yc)^2 < DAC_RADIUS^2
return true
end
end
return false
end

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(is_in_dac))

# ----------------------------------------------------------------------------------
# PHYSICS & FORCING
# ----------------------------------------------------------------------------------

u_vel(x, y, z, t) = 1.0    # m/s wind to the right
v_vel(x, y, z, t) = 0.0
w_vel(x, y, z, t) = 0.0
velocities = (u=u_vel, v=v_vel, w=w_vel)

const CO2_BACKGROUND = 400.0       # ppm
const CO2_PLUME = 420.0
const CO2_TARGET = 400.0
const GAIN = 0.1
const MAX_REMOVAL = 2.0 / 60       # 0.033333... ppm/s
const ENERGY_PER_PPM = 5.0         # units per ppm

@inline function dac_co2_removal(i, j, k, grid, clock, fields)
# Skip solid cells
if grid.immersed_boundary.mask[i, j, k] == 0
return 0.0
end

# Check if near boundary
is_near_boundary = false
for di in -1:1, dj in -1:1
ni, nj = i + di, j + dj
if 1 <= ni <= size(grid, 1) && 1 <= nj <= size(grid, 2)
if grid.immersed_boundary.mask[ni, nj, k] == 0
is_near_boundary = true
break
end
end
end

if is_near_boundary
c = @inbounds fields.c[i, j, k]
if c > CO2_TARGET
sink = -GAIN * (c - CO2_TARGET)
return max(sink, -MAX_REMOVAL)
end
end

return 0.0
end

c_forcing = Forcing(dac_co2_removal, discrete_form=true)

# ----------------------------------------------------------------------------------
# MODEL SETUP
# ----------------------------------------------------------------------------------

model = NonhydrostaticModel(grid=grid,
advection=UpwindBiasedFifthOrder(),
timestepper=:RungeKutta3,
tracers=:c,
velocities=velocities,
forcing=(c=c_forcing,))

set!(model, c = (x, y, z) -> x < 30 ? CO2_PLUME : CO2_BACKGROUND)

# ----------------------------------------------------------------------------------
# SIMULATION SETUP
# ----------------------------------------------------------------------------------

simulation = Simulation(model, Δt=0.5, stop_time=200minutes)

S_field = Field{Center, Center, Center}(model.grid)
global_energy = Float64[]
simulation_times = Float64[]

function energy_callback(sim)
fill!(S_field, 0.0)

grid = sim.model.grid
clock = sim.model.clock
fields = sim.model.tracers

for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
@inbounds S_field[i, j, k] = dac_co2_removal(i, j, k, grid, clock, fields)
end

cell_volume = (Lx / nx) * (Ly / ny) * Lz
total_removal_rate = abs(sum(interior(S_field)) * cell_volume)

energy = ENERGY_PER_PPM * total_removal_rate
push!(global_energy, energy)
push!(simulation_times, sim.clock.time)
end

simulation.output_writers[:co2_output] = JLD2OutputWriter(model, model.tracers,
schedule=TimeInterval(1minute),
filename="co2_dac_feedback.jld2",
overwrite_existing=true)

simulation.callbacks[:progress] = Callback(sim ->
@info @sprintf("Time = %s | Mean CO₂ = %.2f ppm | Max CO₂ = %.2f ppm",
prettytime(sim.clock.time),
mean(interior(sim.model.tracers.c)),
maximum(interior(sim.model.tracers.c))),
TimeInterval(10minutes))

simulation.callbacks[:energy_log] = Callback(energy_callback, IterationInterval(1))

# ----------------------------------------------------------------------------------
# RUN SIMULATION
# ----------------------------------------------------------------------------------

@info "Simulation starting..."
run!(simulation)
@info "Simulation complete."

# ----------------------------------------------------------------------------------
# POST-PROCESSING & VISUALIZATION
# ----------------------------------------------------------------------------------

@info "Generating visualizations..."

filepath = "co2_dac_feedback.jld2"
c_timeseries = FieldTimeSeries(filepath, "c")
times = c_timeseries.times
xc, yc, zc = nodes(c_timeseries[1])
final_time_index = length(times)
c_final = c_timeseries[final_time_index]

function circle_shape(h, k, r)
θ = range(0, 2π, length=100)
return h .+ r * sin.(θ), k .+ r * cos.(θ)
end

# Plot 1: Final concentration
p1 = heatmap(xc, yc, interior(c_final, :, :, 1)',
aspect_ratio=:equal,
xlabel="x (meters)",
ylabel="y (meters)",
title="Final CO₂ Concentration at t = $(prettytime(times[end]))",
c=:thermal,
clims=(CO2_BACKGROUND - 2, CO2_PLUME + 1))

for (xc_dac, yc_dac) in DAC_LOCATIONS
xs, ys = circle_shape(xc_dac, yc_dac, DAC_RADIUS)
plot!(p1, xs, ys, seriestype=:shape, lw=2, c=:white, linecolor=:black,
legend=false, fillalpha=0.2)
end
savefig(p1, "final_co2_concentration.png")

# Plot 2: Energy rate
p2 = plot(simulation_times, global_energy,
xlabel="Time (seconds)",
ylabel="Energy Rate (units/s)",
title="DAC Energy Consumption Rate Over Time",
legend=false,
linewidth=2)
savefig(p2, "dac_energy_consumption.png")

# Animation
anim = @animate for i in 1:length(times)
@info "Animating frame $i of $(length(times))..."
c_t = c_timeseries[i]

p = heatmap(xc, yc, interior(c_t, :, :, 1)',
aspect_ratio=:equal,
xlabel="x (meters)",
ylabel="y (meters)",
title="CO₂ Concentration at t = $(prettytime(times[i])) (ppm)",
c=:thermal,
clims=(CO2_BACKGROUND - 2, CO2_PLUME + 1))

for (xc_dac, yc_dac) in DAC_LOCATIONS
xs, ys = circle_shape(xc_dac, yc_dac, DAC_RADIUS)
plot!(p, xs, ys, seriestype=:shape, lw=2, c=:white, linecolor=:black,
legend=false, fillalpha=0.2)
end

p
end

gif(anim, "co2_dac_evolution.gif", fps=15)

@info "Visualization complete. Files saved: final_co2_concentration.png, dac_energy_consumption.png, co2_dac_evolution.gif"

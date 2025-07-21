using Oceananigans
using Oceananigans.Units: minutes
using Statistics
using Printf
using Plots

# ----------------------
# 1. Model grid setup
# ----------------------
nx, ny = 100, 50        
Lx, Ly = 100.0, 50.0

grid = RectilinearGrid(size=(nx, ny, 1),
                       extent=(Lx, Ly, 1.0),
                       topology=(Bounded, Bounded, Flat))

# ----------------------
# 2. Velocity field (constant rightward wind)
# ----------------------
u(x, y, z, t) = 1.0
v(x, y, z, t) = 0.0
w(x, y, z, t) = 0.0

velocity = (u=u, v=v, w=w)

# ----------------------
# 3. Initialize Oceananigans model
# ----------------------
model = NonhydrostaticModel(grid=grid,
                            advection=UpwindBiasedFifthOrder(),
                            timestepper=:RungeKutta3,
                            tracers=:c,
                            velocities=velocity)

# ----------------------
# 4. Initial CO₂ field: elevated left side plume
# ----------------------
set!(model, c = (x, y, z) -> x < 20 ? 420.0 : 400.0)

# ----------------------
# 5. DAC parameters & locations
# ----------------------
const DAC_LOCATIONS = [(50.0, 25.0), (70.0, 10.0)]
const DAC_RADIUS = 2.0          # meters
const CO2_TARGET = 400.0        # ppm target threshold
const GAIN = 0.2                # proportional gain of removal
const MAX_REMOVAL = 2.0         # ppm per minute max removal rate
const ENERGY_PER_PPM = 5.0      # arbitrary energy units per ppm removed

# ----------------------
# 6. Feedback-controlled DAC sink function
# ----------------------
model.sources.c = ComputedField((x, y, z, t, c) -> begin
    for (xc, yc) in DAC_LOCATIONS
        if abs(x - xc) ≤ DAC_RADIUS && abs(y - yc) ≤ DAC_RADIUS
            if c > CO2_TARGET
                sink = -GAIN * (c - CO2_TARGET)
                return max(sink, -MAX_REMOVAL)  # saturate max removal
            end
        end
    end
    return 0.0
end, dependencies=(model.tracers.c,))

# ----------------------
# 7. Energy tracking container
# ----------------------
global_energy = Float64[]

function energy_callback(sim)
    model = sim.model
    c = model.tracers.c

    energy = 0.0
    for (xc, yc) in DAC_LOCATIONS
        # Find nearest grid index for DAC location
        i = findmin(abs.(grid.xC .- xc))[2]
        j = findmin(abs.(grid.yC .- yc))[2]

        local_c = c[i, j, 1]
        if local_c > CO2_TARGET
            removal = min(GAIN * (local_c - CO2_TARGET), MAX_REMOVAL)
            energy += ENERGY_PER_PPM * removal
        end
    end
    push!(global_energy, energy)
    @info @sprintf("Time = %.2f min | DAC energy use = %.2f units", sim.clock.time / 60, energy)
end

# ----------------------
# 8. Setup and run simulation
# ----------------------
simulation = Simulation(model, Δt=1.0, stop_time=100minutes)

# Output CO₂ every 10 min to disk
simulation.output_writers[:co2_output] = JLD2OutputWriter(model, fields=(:c,),
    schedule=TimeInterval(10minutes), prefix="co2_feedback", force=true)

# Print mean CO₂ every 10 min
simulation.callbacks[:progress] = Callback(every=10) do sim
    cmean = mean(sim.model.tracers.c)
    @info @sprintf("Time = %.2f min | Mean CO₂ = %.2f ppm", sim.clock.time / 60, cmean)
end

# Track energy use every step
simulation.callbacks[:energy_log] = Callback(every=1, func=energy_callback)

# Run simulation
run!(simulation)

# ----------------------
# 9. Visualization (final time)
# ----------------------
using JLD2

file = jldopen("co2_feedback.jld2")
c_final = file["timeseries/c/100.0"]  # CO₂ at 100 min
close(file)

heatmap(grid.xC, grid.yC, c_final[:, :, 1]', aspect_ratio=:equal,
    xlabel="x (m)", ylabel="y (m)", title="CO₂ Concentration at 100 minutes (ppm)",
    c=:thermal, clims=(380, 425))

# Plot DAC energy consumption over time
plot(global_energy, xlabel="Timestep", ylabel="Energy (units)", title="DAC Energy Consumption")

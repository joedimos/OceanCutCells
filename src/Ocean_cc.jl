using Oceananigans, Oceananigans.Units, Printf
using Oceananigans: Grid, VerticallyStretchedRectilinearGrid, architecture
using Oceananigans.Grids: xnodes, ynodes, znodes
using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel
using Oceananigans.TurbulenceClosures: HorizontalLeapfrogDiffusivity, ScalarDiffusivity, nothing_closure
using Oceananigans.Advection: CenteredSecondOrder, UpwindBiasedThirdOrder, nothing_advection
# using Oceananigans.Operators: Δx, ΔzC, ΔzF, KernelFunctionOperation, ConditionalOperation # Keep these operators for reference
using Oceananigans.Fields: ZeroField, Face, Center
using Oceananigans.Forcings: CustomForcing
using Oceananigans.BuoyancyModels: buoyancy_perturbation, SeawaterBuoyancy, LinearEquationOfState
using Oceananigans.Diagnostics: average, create_diagnostic_field
using Oceananigans.BoundaryConditions: fill_halo_regions, FluxBoundaryCondition # Needed for access to halo values

# --- Parameters ---
nx = 40
nz = 20
total_width = 40000.0
total_depth = 3000.0
x_domain = (0, total_width)
z_domain = (-total_depth, 0)

dt = 100.0           # Time step
Kh = 5.0; Kv = 1e-4    # Tracer diffusion
Ah = 20.0; Av = 1e-3   # Momentum viscosity (Used in vertical diffusion of u)
rho0 = 1025.0; g = 9.81
alpha = 2e-4; beta = 7.4e-4; T0 = 15.0; S0 = 35.5
wind_stress_tau_x = 0.05 # Apply eastward wind stress (N/m^2)
bottom_drag_coeff = 2e-3  # Bottom drag coeff

# Define a bathymetry function and its derivative
function bathymetry_profile(x)
    base = 500 + 2000 * (x / total_width)
    xc = total_width * 0.6
    w = total_width * 0.08
    ridge = 1500 * exp(-((x - xc)^2) / (2 * w^2))
    depth = base - ridge
    return max(50.0, depth) # Min depth
end

function bathymetry_derivative(x)
    dbase_dx = 2000.0 / total_width
    xc = total_width * 0.6
    w = total_width * 0.08
    # Derivative of exp(-((x - xc)^2) / (2 * w^2)) with respect to x is
    # exp(-((x - xc)^2) / (2 * w^2)) * (-(x - xc) / w^2)
    dridge_dx = 1500.0 * exp(-((x - xc)^2) / (2 * w^2)) * (-(x - xc) / w^2)
    return dbase_dx - dridge_dx
end


# --- Grid Generation ---
# Non-uniform z-spacing: More resolution near the surface
z_stretching = 1.3 # Factor for vertical stretching
z_faces_func(k) = -total_depth * (k/nz)^z_stretching

# Create a VerticallyStretchedRectilinearGrid
grid = VerticallyStretchedRectilinearGrid(
    topology = (Bounded, Flat, Bounded),
    size = (nx, nz),
    x = x_domain,
    z_faces = z_faces_func.(0:nz),
    halo = (2, 2) # 2 in x, 2 in z
)

# --- Cut Cell Implementation: hFac fields and Wet Areas/Volumes (Initial) ---

# Define hFacC, hFacW, and hFacS fields
hFacC_initial = Field{Center, Center, Center}(grid)
hFacW_initial = Field{Face, Center, Center}(grid)
hFacS_initial = Field{Center, Center, Face}(grid)

# Calculate initial hFac values based on bathymetry
x_centers = xnodes(hFacC_initial, Center(), Center(), Center())[:, 1, 1]
bathymetry_c = bathymetry_profile.(x_centers) # At C locations

for i in 1:grid.Nx, k in 1:grid.Nz
    z_top_face = grid.zᵃᵃᶠ[k]
    z_bot_face = grid.zᵃᵃᶠ[k+1]
    dz_k = grid.Δzᵃᵃᶜ[k]

    bathy_z_at_c = -bathymetry_c[i]

    if z_top_face < bathy_z_at_c
        hFacC_initial.data[i, 1, k] = 0.0
    elseif z_bot_face >= bathy_z_at_c
        hFacC_initial.data[i, 1, k] = 1.0
    else # Partial cell
        hFacC_initial.data[i, 1, k] = clamp((bathy_z_at_c - z_bot_face) / dz_k, 0.0, 1.0)
    end
end

# Simple connectivity masks for faces based on initial hFacC
for i in 1:grid.Nx+1, k in 1:grid.Nz
    if i > grid.Nx || i == 1 # Solid West/East boundaries
         hFacW_initial.data[i, 1, k] = 0.0
         continue
    end
    is_wet_im1k = hFacC_initial.data[i-1, 1, k] > 1e-10
    is_wet_ik = hFacC_initial.data[i, 1, k] > 1e-10
    hFacW_initial.data[i, 1, k] = (is_wet_im1k && is_wet_ik) ? 1.0 : 0.0
end

for i in 1:grid.Nx, k in 1:grid.Nz+1
     if k > grid.Nz # Solid bottom boundary
        hFacS_initial.data[i, 1, k] = 0.0
        continue
     end
     if k == 1 # Open surface (z=0)
         hFacS_initial.data[i, 1, k] = hFacC_initial.data[i, 1, k] > 1e-10 ? 1.0 : 0.0
         continue
     end
    is_wet_ikm1 = hFacC_initial.data[i, 1, k-1] > 1e-10
    is_wet_ik = hFacC_initial.data[i, 1, k] > 1e-10
    hFacS_initial.data[i, 1, k] = (is_wet_ikm1 && is_wet_ik) ? 1.0 : 0.0
end

# Apply hard cut-off
hFacC_initial[hFacC_initial.data .< 1e-10] .= 0.0
hFacW_initial[hFacW_initial.data .< 1e-10] .= 0.0
hFacS_initial[hFacS_initial.data .< 1e-10] .= 0.0


# --- Calculate Initial Wet Areas and Volumes (before merging) ---

wet_cell_volume_initial = Field{Center, Center, Center}(grid)
wet_face_area_x_initial = Field{Face, Center, Center}(grid)
wet_face_area_z_initial = Field{Center, Center, Face}(grid)

# Wet volume (at C)
@inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
    wet_cell_volume_initial.data[i+grid.Hx, j+grid.Hy, k+grid.Hz] = hFacC_initial[i, j, k] * grid.Δxᶜᵃᵃ[i] * grid.Δyᵃᶜᵃ[j] * grid.Δzᵃᵃᶜ[k]
end

# Wet face area x (at F, west face i)
@inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
    wet_face_area_x_initial.data[i+grid.Hx-1, j+grid.Hy, k+grid.Hz] = hFacW_initial[i, j, k] * grid.Δyᵃᶜᵃ[j] * grid.Δzᵃᵃᶜ[k]
end

# Wet face area z (at F, bottom face k)
@inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz+1
     wet_face_area_z_initial.data[i+grid.Hx, j+grid.Hy, k+grid.Hz-1] = hFacS_initial[i, j, k] * grid.Δxᶜᵃᵃ[i] * grid.Δyᵃᶜᵃ[j]
end


# --- Implement Cell Merging ---

@info "Performing cell merging..."

# Copies to modify
hFacC_merged = deepcopy(hFacC_initial)
wet_cell_volume_merged = deepcopy(wet_cell_volume_initial)
# NOT merging face areas for now (wet_face_area_x_merged, wet_face_area_z_merged will be deepcopies initially)
wet_face_area_x_merged = deepcopy(wet_face_area_x_initial)
wet_face_area_z_merged = deepcopy(wet_face_area_z_initial)

# Bathymetry derivative at cell centers in x
dhdx_c = bathymetry_derivative.(x_centers)

# Mask for cells that will be merged away (source cells)
source_mask = fill(false, grid.Nx, grid.Ny, grid.Nz)
# Store the target index for each source cell
target_idx = fill((-1, -1, -1), grid.Nx, grid.Ny, grid.Nz) # Store (i,j,k) tuple

# Identify source cells and their targets
@inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
    regular_cell_volume = grid.Δxᶜᵃᵃ[i] * grid.Δyᵃᶜᵃ[j] * grid.Δzᵃᵃᶜ[k]
    is_small_cut_cell = (hFacC_initial[i, j, k] > 1e-10) && (wet_cell_volume_initial[i, j, k] < 0.5 * regular_cell_volume)

    if is_small_cut_cell
        # Determine merge direction based on dh/dx (2D logic)
        grad_x = dhdx_c[i]

        i_t, j_t, k_t = i, j, k # Default target is self

        if abs(grad_x) <= 1.0
            # Vertical merge upwards
            k_t = k - 1
        elseif grad_x > 1.0
            # Horizontal merge left (-x)
            i_t = i - 1
        elseif grad_x < -1.0
            # Horizontal merge right (+x)
            i_t = i + 1
        end

        # Check if the target cell is valid and wet (within interior domain)
        is_target_valid = (i_t >= 1 && i_t <= grid.Nx && j_t >= 1 && j_t <= grid.Ny && k_t >= 1 && k_t <= grid.Nz)
        is_target_wet = is_target_valid && (hFacC_initial[i_t, j_t, k_t] > 1e-10)

        if is_target_wet
            # Mark this cell as a source cell
            source_mask[i, j, k] = true
            target_idx[i, j, k] = (i_t, j_t, k_t)
        else
             @warn "Small cell at ($i, $j, $k) could not be merged (invalid or dry target)."
             # This cell remains a small cell. It will be computed on but might cause issues.
             # In a real implementation, more sophisticated handling or error reporting needed.
             # For now, we just leave it as is if it can't merge.
             source_mask[i, j, k] = false # Do not merge this cell
        end
    end
end

# Perform the merging of volumes
@inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
    if source_mask[i, j, k]
        (i_t, j_t, k_t) = target_idx[i, j, k]

        # Add source volume to target volume
        wet_cell_volume_merged.data[i_t+grid.Hx-1, j_t+grid.Hy-1, k_t+grid.Hz-1] += wet_cell_volume_initial[i, j, k]

        # Set source volume to zero
        wet_cell_volume_merged.data[i+grid.Hx-1, j+grid.Hy-1, k+grid.Hz-1] = 0.0

        # Mark source hFacC as zero
        hFacC_merged.data[i+grid.Hx-1, j+grid.Hy-1, k+grid.Hz-1] = 0.0

        # Note: Face areas are NOT merged here. This is a simplification.
    end
end

@info "Cell merging complete. Number of merged cells: $(sum(source_mask))"


# --- Calculate Reciprocal Wet Volume for Merged Grid ---
recip_wet_cell_volume_merged = Field{Center, Center, Center}(grid)
@inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
    if wet_cell_volume_merged[i,j,k] > 1e-15 # Use the modified volume
        recip_wet_cell_volume_merged.data[i+grid.Hx, j+grid.Hy, k+grid.Hz] = 1.0 / wet_cell_volume_merged[i,j,k]
    else
        recip_wet_cell_volume_merged.data[i+grid.Hx, j+grid.Hy, k+grid.Hz] = 0.0 # Reciprocal is 0 for dry cells (merged-away cells)
    end
end

# Ensure halo regions are filled for merged fields
fill_halo_regions!(hFacC_merged)
fill_halo_regions!(wet_cell_volume_merged)
fill_halo_regions!(recip_wet_cell_volume_merged)
fill_halo_regions!(wet_face_area_x_merged) # These are still initial areas
fill_halo_regions!(wet_face_area_z_merged) # These are still initial areas


# --- Bundle Cut-Cell Parameters (using merged geometry) ---
struct CutCellParameters
    hFacC :: Field{Center, Center, Center}
    hFacW :: Field{Face, Center, Center}
    hFacS :: Field{Center, Center, Face}
    wet_cell_volume   :: Field{Center, Center, Center}
    recip_wet_cell_volume :: Field{Center, Center, Center}
    wet_face_area_x :: Field{Face, Center, Center} # At West face i
    wet_face_area_z :: Field{Center, Center, Face} # At Bottom face k
    # Add coefficients
    rho0 :: Float64
    g :: Float64
    Kh :: Float64
    Kv :: Float64
    Ah :: Float64
    Av :: Float64
    bottom_drag_coeff :: Float64
end

# Use the _merged fields in the parameters
cut_cell_params = CutCellParameters(hFacC_merged, hFacW_initial, hFacS_initial, # hFacW/S are not merged in this approach
                                     wet_cell_volume_merged, recip_wet_cell_volume_merged,
                                     wet_face_area_x_merged, wet_face_area_z_merged, # Wet face areas are not merged
                                     rho0, g, Kh, Kv, Ah, Av, bottom_drag_coeff)


# --- Cut-Cell Aware Forcing Functions (Flux Based) ---
# These remain largely the same, but use the merged geometric fields
# from cut_cell_params. The key change is that `wet_cell_volume` and `recip_wet_cell_volume`
# will have the merged values, and `hFacC` will be zero for merged-away cells.
# The forcings use `hFacC` to gate the tendency computation (`if cc_params.hFacC[i, j, k] > 1e-10`).

# 1. Pressure Gradient Force for u (at Face, Center, Center)
# This forcing accesses p (C), hFacW (F), hFacC (C), wet_face_area_x (F), recip_wet_cell_volume (C)?
# PGF adds force to u-momentum (F). Tendency is F/m.
# This forcing is added to `∂u∂t` which is acceleration (Force/Mass).
# PGF = -(1/rho0) * ∇p * V_u_cell
# div(stress)/rho0 = (1/rho0) * (stress_above - stress_below)/dz_u_cell
# Advection adds -(u.∇)u * V_u_cell
# Oh, I made a mistake in the previous forcings. They should be adding ACCELERATION (Force/Mass), not Force.
# PGF: -(1/rho0) * ∂p/∂x.  The flux divergence is d/dx(p). The discrete form is (p_i - p_im1)/dx.
# This discrete gradient is the acceleration.

# Let's correct the forcing functions to add acceleration directly where possible,
# or recalculate flux divergence then divide by effective volume if necessary for cut cells.

# Re-evaluating the custom forcings based on Oceananigans' approach to CustomForcing:
# CustomForcing functions add to the TENDENCY field.
# For tracer T (at C), tendency is ∂T/∂t (Units T/s). Flux divergence needs to be divided by Wet Cell Volume.
# For momentum u (at F), tendency is ∂u/∂t (Units m/s²). Forces (PGF, Drag, etc.) are already m/s². Advection (u.∇)u is m/s². Diffusion Av * ∇²u is m/s².

# Let's rewrite the forcings more carefully.

# 1. Pressure Gradient Force for u (at Face, Center, Center)
# Adds ∂u∂t = -(1/ρ₀) ∂p/∂x
# Uses Oceananigans' built-in gradient operator for discrete (p_i - p_im1)/dx at F,C,C.
# Need to mask this gradient by hFacW.
function add_cut_cell_pressure_gradient_force!(∂u∂t, model)
    grid = model.grid
    p = model.pressure
    cc_params = model.parameters.cut_cell_params

    # Use Oceananigans' discrete gradient operator ∂xᶠᶜᶜ
    # This operator is defined over the internal Face points (1:Nx+1 in x)
    # It computes (p[i,j,k] - p[i-1,j,k]) / grid.Δxᶠᶜᶜ[i,j,k]
    p_gradient = ∂xᶠᶜᶜ(1.0 * p) # Compute gradient as a field

    # Add -(1/rho0) * p_gradient to ∂u∂t, masked by hFacW
    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
        # hFacW gates the flux across the face. If the face is solid (hFacW=0), the PGF should be balanced by the wall.
        # If the face is wet, PGF contributes to acceleration.
        if cc_params.hFacW[i, j, k] > 1e-10
             # Masking by hFacW is crucial for the force term.
             # The PGF itself is applied to the u-point. The hFacW indicates if the u-point is "active".
             ∂u∂t.data[i+grid.Hx-1, j+grid.Hy, k+grid.Hz] -= (1.0 / cc_params.rho0) * p_gradient[i, j, k]
        else
             # If hFacW is zero, the face is solid. The tendency should be zero here implicitly due to BCs on u.
             # But explicitly setting the forcing to zero ensures it doesn't contribute if hFacW is exactly zero.
             ∂u∂t.data[i+grid.Hx-1, j+grid.Hy, k+grid.Hz] += 0.0 # Explicitly add zero for clarity
        end
    end
    return nothing
end


# 2. Tracer Advection for T (at Center, Center, Center) - Flux-Based Upwind
# Adds ∂T/∂t = - (1/V_cell) * ∇ . (u T) * V_cell = - ∇ . (u T)
# Uses CustomForcing to calculate flux divergence and divide by Wet Cell Volume.
function add_cut_cell_advection_T!(∂T∂t, model)
    grid = model.grid
    T = model.tracers.T
    # w is the diagnosed w field, must be passed as a field or accessed from parameters
    w_field = model.velocities.w # Access the velocity field updated by diagnostic

    cc_params = model.parameters.cut_cell_params

    # Loop over the interior grid points where T tendency is computed (C,C,C)
    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
        # Only compute tendency for wet cells (using the merged hFacC mask)
        if cc_params.hFacC[i, j, k] > 1e-10

            # --- Horizontal Fluxes (East - West) of T ---
            # Flux across west face (i) - where u[i,j,k] lives. Check hFacW[i,j,k].
            # Flux = u_face * T_upstream * WetArea_face
            flux_w_T_horiz = 0.0
            if cc_params.hFacW[i, j, k] > 1e-10
                 # T_upstream: If u[i,j,k] > 0, source is cell (i-1,k). If u[i,j,k] <= 0, source is cell (i,k).
                 # Need to check if source cells are wet.
                 is_wet_im1k = (i > 1 && cc_params.hFacC[i-1, j, k] > 1e-10)
                 is_wet_ik = (i <= grid.Nx && cc_params.hFacC[i, j, k] > 1e-10) # Current cell is wet

                 if is_wet_im1k || is_wet_ik # Flux calculated if face is wet and connects to at least one wet cell
                     u_face = u[i, j, k]
                     if u_face > 0 # Upstream is (i-1, j, k)
                         T_upstream = is_wet_im1k ? T[i-1, j, k] : 0.0 # If upstream cell is dry, use 0
                     else # Upstream is (i, j, k)
                         T_upstream = is_wet_ik ? T[i, j, k] : 0.0 # If upstream cell is dry, use 0 (should be wet here)
                     end
                     flux_w_T_horiz = u_face * T_upstream * cc_params.wet_face_area_x[i, j, k]
                 end
            end

            # Flux across east face (i+1) - where u[i+1,j,k] lives. Check hFacW[i+1,j,k].
            flux_e_T_horiz = 0.0
            if (i < grid.Nx) && cc_params.hFacW[i+1, j, k] > 1e-10
                 is_wet_ik = (i <= grid.Nx && cc_params.hFacC[i, j, k] > 1e-10) # Current cell is wet
                 is_wet_ip1k = (i < grid.Nx && cc_params.hFacC[i+1, j, k] > 1e-10)

                 if is_wet_ik || is_wet_ip1k
                    u_face = u[i+1, j, k]
                    if u_face > 0 # Upstream is (i, j, k)
                        T_upstream = is_wet_ik ? T[i, j, k] : 0.0 # If upstream cell is dry, use 0 (should be wet)
                    else # Upstream is (i+1, j, k)
                        T_upstream = is_wet_ip1k ? T[i+1, j, k] : 0.0
                    end
                    flux_e_T_horiz = u_face * T_upstream * cc_params.wet_face_area_x[i+1, j, k]
                 end
            end

            # Horizontal Divergence = (East Flux - West Flux)
            horiz_flux_div_T = flux_e_T_horiz - flux_w_T_horiz

            # --- Vertical Fluxs (North - South) of T ---
            # Flux across south face (k) - where w[i,j,k] lives. Check hFacS[i,j,k].
            # Flux = w_face * T_upstream * WetArea_face
            flux_s_T_vert = 0.0
            if cc_params.hFacS[i, j, k] > 1e-10
                 # T_upstream: If w[i,j,k] > 0, source is cell (i,k+1). If w[i,j,k] <= 0, source is cell (i,k).
                 is_wet_ik = (k <= grid.Nz && cc_params.hFacC[i, j, k] > 1e-10) # Current cell is wet
                 is_wet_ikp1 = (k < grid.Nz && cc_params.hFacC[i, j, k+1] > 1e-10)

                 if is_wet_ik || is_wet_ikp1
                     w_face = w_field[i, j, k]
                     if w_face > 0 # Upstream is (i, j, k+1)
                         T_upstream = is_wet_ikp1 ? T[i, j, k+1] : 0.0
                     else # Upstream is (i, j, k)
                         T_upstream = is_wet_ik ? T[i, j, k] : 0.0 # Should be wet
                     end
                     flux_s_T_vert = w_face * T_upstream * cc_params.wet_face_area_z[i, j, k]
                 end
            end

            # Flux across north face (k-1) - where w[i,j,k-1] lives. Check hFacS[i,j,k-1].
            flux_n_T_vert = 0.0
            if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10
                  is_wet_ikm1 = (k > 1 && cc_params.hFacC[i, j, k-1] > 1e-10)
                  is_wet_ik = (k <= grid.Nz && cc_params.hFacC[i, j, k] > 1e-10) # Current cell is wet

                  if is_wet_ikm1 || is_wet_ik
                     w_face = w_field[i, j, k-1]
                     if w_face > 0 # Upstream is (i, j, k)
                         T_upstream = is_wet_ik ? T[i, j, k] : 0.0 # Should be wet
                     else # Upstream is (i, j, k-1)
                         T_upstream = is_wet_ikm1 ? T[i, j, k-1] : 0.0
                     end
                     flux_n_T_vert = w_face * T_upstream * cc_params.wet_face_area_z[i, j, k-1]
                  end
            end

            # Vertical Divergence = (North Flux - South Flux)
            vert_flux_div_T = flux_n_T_vert - flux_s_T_vert

            # Total Divergence
            total_flux_div_T = horiz_flux_div_T + vert_flux_div_T

            # Tendency = - Total Divergence / Wet Cell Volume (Use merged volume)
            ∂T∂t[i, j, k] -= total_flux_div_T * cc_params.recip_wet_cell_volume[i, j, k]

        end # if wet cell (using merged hFacC)
    end # Loop
    return nothing
end


# 3. Tracer Diffusion for T (at Center, Center, Center) - Flux-Based
# Adds ∂T/∂t = - (1/V_cell) * ∇ . (-K ∇ T) * V_cell = ∇ . (K ∇ T)
function add_cut_cell_diffusion_T!(∂T∂t, model)
    grid = model.grid
    T = model.tracers.T
    cc_params = model.parameters.cut_cell_params

    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
        # Only compute tendency for wet cells (using merged hFacC)
        if cc_params.hFacC[i, j, k] > 1e-10

            # --- Horizontal Fluxes (East - West) of T Gradient ---
            # Flux across west face (i) = -Kh * grad_x(T) * WetArea_face_x
            flux_w_T_horiz = 0.0
            if cc_params.hFacW[i, j, k] > 1e-10 # Face must be wet
                 # Gradient requires both adjacent cells to be wet (using merged hFacC)
                 is_wet_im1k = (i > 1 && cc_params.hFacC[i-1, j, k] > 1e-10)
                 is_wet_ik = (i <= grid.Nx && cc_params.hFacC[i, j, k] > 1e-10) # Current cell is wet

                 if is_wet_im1k && is_wet_ik # Gradient calculated only if face connects two wet cells
                     Δx_centers = grid.Δxᶠᶜᶜ[i, j, k] # Distance between T at (i-1,k) and (i,k)
                      if Δx_centers > 1e-10
                         grad = (T[i, j, k] - T[i-1, j, k]) / Δx_centers
                         flux_w_T_horiz = -cc_params.Kh * grad * cc_params.wet_face_area_x[i, j, k]
                      end
                 end
            end

            # Flux across east face (i+1) = -Kh * grad_x(T) * WetArea_face_x
            flux_e_T_horiz = 0.0
            if (i < grid.Nx) && cc_params.hFacW[i+1, j, k] > 1e-10
                 is_wet_ik = (i <= grid.Nx && cc_params.hFacC[i, j, k] > 1e-10) # Current cell is wet
                 is_wet_ip1k = (i < grid.Nx && cc_params.hFacC[i+1, j, k] > 1e-10)

                 if is_wet_ik && is_wet_ip1k
                    Δx_centers = grid.Δxᶠᶜᶜ[i+1, j, k]
                    if Δx_centers > 1e-10
                       grad = (T[i+1, j, k] - T[i, j, k]) / Δx_centers
                       flux_e_T_horiz = -cc_params.Kh * grad * cc_params.wet_face_area_x[i+1, j, k]
                    end
                 end
            end

            # Horizontal Flux Divergence = (East Flux - West Flux)
            horiz_flux_div_T = flux_e_T_horiz - flux_w_T_horiz

            # --- Vertical Fluxes (North - South) of T Gradient ---
            # Flux across south face (k) = -Kv * grad_z(T) * WetArea_face_z
            flux_s_T_vert = 0.0
            if cc_params.hFacS[i, j, k] > 1e-10 # Face must be wet
                 # Gradient requires both adjacent cells to be wet (using merged hFacC)
                 is_wet_ik = (k <= grid.Nz && cc_params.hFacC[i, j, k] > 1e-10) # Current cell is wet
                 is_wet_ikp1 = (k < grid.Nz && cc_params.hFacC[i, j, k+1] > 1e-10)

                 if is_wet_ik && is_wet_ikp1
                     Δz_centers = grid.Δzᵃᵃᶠ[k, j, i] # Distance between T at (i,k+1) and (i,k)
                      if Δz_centers > 1e-10
                         grad = (T[i, j, k+1] - T[i, j, k]) / Δz_centers
                         # Flux is positive *upwards* if Kv*grad is positive (T increases downwards -> heat flux upwards)
                         # Flux *across face k* (which is the SOUTH face of cell k) is directed downwards
                         flux_s_T_vert = -cc_params.Kv * grad * cc_params.wet_face_area_z[i, j, k]
                      end
                 end
            end

            # Flux across north face (k-1) = -Kv * grad_z(T) * WetArea_face_z
            flux_n_T_vert = 0.0
            if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10 # Face must be wet
                  is_wet_ikm1 = (k > 1 && cc_params.hFacC[i, j, k-1] > 1e-10)
                  is_wet_ik = (k <= grid.Nz && cc_params.hFacC[i, j, k] > 1e-10) # Current cell is wet

                  if is_wet_ikm1 && is_wet_ik
                     Δz_centers = grid.Δzᵃᵃᶠ[k-1, j, i] # Distance between T at (i,k-1) and (i,k)
                     if Δz_centers > 1e-10
                         grad = (T[i, j, k] - T[i, j, k-1]) / Δz_centers # Gradient between k-1 and k
                         # Flux across face k-1 (NORTH face of cell k) is directed upwards
                         flux_n_T_vert = -cc_params.Kv * grad * cc_params.wet_face_area_z[i, j, k-1]
                     end
                  end
            end

            # Vertical Flux Divergence = (North Flux - South Flux)
            vert_flux_div_T = flux_n_T_vert - flux_s_T_vert

            # Total Flux Divergence
            total_flux_div_T = horiz_flux_div_T + vert_flux_div_T

            # Tendency = - Total Divergence / Wet Cell Volume (Use merged volume)
            ∂T∂t[i, j, k] -= total_flux_div_T * cc_params.recip_wet_cell_volume[i, j, k]

        end # if wet cell (using merged hFacC)
    end # Loop
    return nothing
end


# 4. Vertical Advection for u (at Face, Center, Center) - Flux-Based Upwind
# Adds ∂u/∂t = - (1/V_u_cell) * ∇ . (w u) * V_u_cell = - ∇ . (w u)
# Volume of u cell is approximately Δx * Δy * Δzᵃᵃᶜ[k] where u[i,j,k] lives.
# The divergence term for vertical advection ∂/∂z (w u) contributes to ∂u/∂t.
# This term is calculated at (F,C,C). The flux w*u needs to be evaluated at (F,C,F) locations.
# w is at (C,C,F). u is at (F,C,C).
# w needs to be interpolated to (F,C,F). u needs to be interpolated to (C,C,F).
# This is complex. Let's use a simplified approach: Sample w at the face center (C,C,F) and u by upwinding from (F,C,C).
# The area is wet_face_area_z (at C,C,F). Let's use this as an approximation.
# The divergence is calculated by (Flux_top - Flux_bottom) / Δz_u_cell.
function add_cut_cell_vertical_advection_u!(∂u∂t, model)
    grid = model.grid
    u = model.velocities.u
    w_field = model.velocities.w # w field updated by diagnostic
    cc_params = model.parameters.cut_cell_params

    # Loop over the interior domain where u-tendency is calculated (F,C,C)
    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
         # Only compute tendency for wet u-faces (using original hFacW)
         if cc_params.hFacW[i, j, k] > 1e-10

            # --- Vertical Fluxes (North - South) of u ---
            # Flux across horizontal face at z_F[k+1] (below u[i,j,k]) - where w[i,j,k] lives. Check hFacS[i,j,k].
            # This flux is w * u_face * WetArea_face_z
            # w is at (C,C,F), u is at (F,C,C), WetArea_face_z is at (C,C,F)
            # Approximation: use w[i,j,k] and wet_face_area_z[i,j,k] for the flux at this face location.
            # Upwind u from u[i,j,k] and u[i,j,k+1].
            flux_below_u_vert = 0.0 # Flux out the bottom of u[i,j,k] cell
            if cc_params.hFacS[i, j, k] > 1e-10 # Face must be wet

                 # Check connectivity of adjacent u-points for upwinding u using original hFacW
                 is_wet_u_ik = (k <= grid.Nz && cc_params.hFacW[i, j, k] > 1e-10) # u[i,j,k] wet? (Should be true if this loop iteration is reached)
                 is_wet_u_ikp1 = (k < grid.Nz && cc_params.hFacW[i, j, k+1] > 1e-10) # u[i,j,k+1] wet?

                 if is_wet_u_ik || is_wet_u_ikp1 # Flux calculated if face is wet and connects to at least one wet u-point
                     w_face = w_field[i, j, k] # w velocity at the face
                     if w_face > 0 # Upstream is u[i,j,k+1]
                         u_upstream = is_wet_u_ikp1 ? u[i, j, k+1] : 0.0 # If upstream u is dry, use 0
                     else # Upstream is u[i,j,k]
                         u_upstream = is_wet_u_ik ? u[i, j, k] : 0.0 # Should be wet
                     end
                     flux_below_u_vert = w_face * u_upstream * cc_params.wet_face_area_z[i, j, k]
                 end
            end

            # Flux across horizontal face at z_F[k] (above u[i,j,k]) - where w[i,j,k-1] lives. Check hFacS[i,j,k-1].
            flux_above_u_vert = 0.0 # Flux into the top of u[i,j,k] cell
            if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10 # Face must be wet

                  is_wet_u_ikm1 = (k > 1 && cc_params.hFacW[i, j, k-1] > 1e-10) # u[i,j,k-1] wet?
                  is_wet_u_ik = (k <= grid.Nz && cc_params.hFacW[i, j, k] > 1e-10) # u[i,j,k] wet? (Should be true)

                   if is_wet_u_ikm1 || is_wet_u_ik # Flux calculated if face is wet and connects to at least one wet u-point
                     w_face = w_field[i, j, k-1]
                     if w_face > 0 # Upstream is u[i,j,k]
                         u_upstream = is_wet_u_ik ? u[i, j, k] : 0.0 # Should be wet
                     else # Upstream is u[i,j,k-1]
                         u_upstream = is_wet_u_ikm1 ? u[i, j, k-1] : 0.0
                     end
                     flux_above_u_vert = w_face * u_upstream * cc_params.wet_face_area_z[i, j, k-1]
                   end
            end

            # Vertical Divergence = (Flux_into_top - Flux_out_bottom) / Volume_height
            # Volume_height is the height of the u-cell control volume, approx Δzᵃᵃᶜ[k].
            Δz_u_cell = grid.Δzᵃᵃᶜ[k]
            if Δz_u_cell > 1e-10
                 # Divergence of vertical flux = (Flux_above - Flux_below) / Δz_u_cell
                 vert_flux_div_u = (flux_above_u_vert - flux_below_u_vert) / Δz_u_cell
                 ∂u∂t.data[i+grid.Hx-1, j+grid.Hy, k+grid.Hz] -= vert_flux_div_u # Subtract divergence for tendency
            end
         end # If hFacW > 0 (using original hFacW)
    end # Loop
    return nothing
end


# 5. Vertical Diffusion for u (at Face, Center, Center) - Flux-Based Vertical Viscosity
# Adds ∂u/∂t = ∂/∂z (Av ∂u/∂z)
# Uses a centered stencil approximation for the Laplacian, masked by hFacW connectivity.
function add_cut_cell_vertical_diffusion_u!(∂u∂t, model)
    grid = model.grid
    u = model.velocities.u
    cc_params = model.parameters.cut_cell_params

    # Loop over the interior domain where u-tendency is calculated (F,C,C)
    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
         # Only compute tendency for wet u-faces (using original hFacW)
         if cc_params.hFacW[i, j, k] > 1e-10

            # Centered difference approximation for ∂/∂z (Av ∂u/∂z) at (i,j,k)
            # Av * [ (u[i,j,k+1] - u[i,j,k]) / Δz_face_k+1 - (u[i,j,k] - u[i,j,k-1]) / Δz_face_k ] / Δz_cell_k
            # Δz_face_k+1 is distance between u[k] and u[k+1] (approx Δzᵃᵃᶠ[k] )
            # Δz_face_k is distance between u[k-1] and u[k] (approx Δzᵃᵃᶠ[k-1])
            # Δz_cell_k is height of u[i,j,k] cell (approx Δzᵃᵃᶜ[k])

            # Need to check connectivity using original hFacW
            is_wet_u_ik = (k <= grid.Nz && cc_params.hFacW[i, j, k] > 1e-10) # u[i,j,k] wet? (Should be true)
            is_wet_u_ikm1 = (k > 1 && cc_params.hFacW[i, j, k-1] > 1e-10) # u[i,j,k-1] wet?
            is_wet_u_ikp1 = (k < grid.Nz && cc_params.hFacW[i, j, k+1] > 1e-10) # u[i,j,k+1] wet?

            # Vertical gradient below u[i,j,k] (between k and k+1)
            grad_below = 0.0
            if is_wet_u_ik && is_wet_u_ikp1 # Need both points to calculate gradient
                Δz_between = grid.Δzᵃᵃᶠ[k, j, i]
                 if Δz_between > 1e-10
                    grad_below = (u[i, j, k+1] - u[i, j, k]) / Δz_between
                 end
            end # If points not connected, gradient is effectively zero

            # Vertical gradient above u[i,j,k] (between k-1 and k)
            grad_above = 0.0
             if is_wet_u_ikm1 && is_wet_u_ik # Need both points
                 Δz_between = grid.Δzᵃᵃᶠ[k-1, j, i]
                 if Δz_between > 1e-10
                     grad_above = (u[i, j, k] - u[i, j, k-1]) / Δz_between
                 end
             end # If points not connected, gradient is effectively zero

            # Vertical Laplacian: Divergence of gradient
            Δz_u_cell = grid.Δzᵃᵃᶜ[k]
            if Δz_u_cell > 1e-10
                # Divergence of flux = (Flux_above - Flux_below) / Volume_height
                # Flux is -Av * grad. So div(-Av grad) = -Av * div(grad) = -Av * Laplacian
                # The term added to ∂u∂t is + Av * Laplacian
                # Laplacian approx = (grad_above - grad_below) / Δz_u_cell
                laplacian_z = (grad_above - grad_below) / Δz_u_cell
                ∂u∂t.data[i+grid.Hx-1, j+grid.Hy, k+grid.Hz] += cc_params.Av * laplacian_z
            end

         end # If hFacW > 0 (using original hFacW)
    end # Loop
    return nothing
end


# 6. Bottom Drag (at Face, Center, Center)
# Adds ∂u/∂t = - Cd |u| u / Δz_bottom_cell
function add_cut_cell_bottom_drag!(∂u∂t, model)
    grid = model.grid
    u = model.velocities.u
    cc_params = model.parameters.cut_cell_params # Access bundled parameters

    Cd = cc_params.bottom_drag_coeff

    # Apply drag to the bottom-most wet u-point in each column
    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny
        k_bottom = -1
        # Find the deepest wet u-face point [i,j,k]. Iterate from bottom up.
        # Use original hFacW to find wet u-points.
        for k in grid.Nz : -1 : 1
           if cc_params.hFacW[i, j, k] > 1e-10
               k_bottom = k
               break # Found the bottom-most wet u point in this column (i,j)
           end
        end

        if k_bottom != -1
            u_bottom = u[i, j, k_bottom]
            # Drag force per unit mass = - Cd * |u_bottom| * u_bottom / H_eff
            # H_eff is the vertical extent over which the drag is distributed. Use the thickness of the cell: Δzᵃᵃᶜ[k_bottom].
            H_eff = grid.Δzᵃᵃᶜ[k_bottom]

            if H_eff > 1e-10
                drag_accel = - Cd * abs(u_bottom) * u_bottom / H_eff
                ∂u∂t.data[i+grid.Hx-1, j+grid.Hy, k_bottom+grid.Hz] += drag_accel
            end
        end
    end
    return nothing
end


# --- Custom Diagnostics ---
# Create a field for the diagnosed w
cut_cell_w = create_diagnostic_field(model.velocities.w)

# Continuity equation diagnostic to calculate w (at C,C,F) from u (at F,C,C)
# ∇ . u = 0
# ∂u/∂x + ∂w/∂z = 0 (in 2D x-z)
# Integrate ∂w/∂z = -∂u/∂x from bottom up to get w.
# w[k] - w[k+1] = - ∫(∂u/∂x) dz approx - (∂u/∂x)_cell * Δz_cell
# w[i,j,k] * WetArea_z[i,j,k] = w[i,j,k+1] * WetArea_z[i,j,k+1] + Net_Horiz_Vol_Flux_INTO_cell(i,j,k)
# Net Horiz Vol Flux INTO cell (i,j,k) = u_West_face * WetArea_x_West - u_East_face * WetArea_x_East
# u_West_face = u[i,j,k], Area_West = wet_face_area_x[i,j,k]
# u_East_face = u[i+1,j,k], Area_East = wet_face_area_x[i+1,j,k]
# This is the same continuity diagnostic as before, just ensuring it uses the (non-merged) face areas.
function diagnose_cut_cell_w!(w_field, model)
    grid = model.grid
    u = model.velocities.u
    cc_params = model.parameters.cut_cell_params

    fill_halo_regions!(u) # Ensure u halos are up-to-date

    # Array to store vertical volume flux (w*Area)
    vertical_volume_flux = deepcopy(w_field)
    vertical_volume_flux.data .= 0.0

    # Integrate from bottom up (k = grid.Nz to 1)
    # The vertical face k of cell (i,j,k) is at z_F[k+1] if z_F is decreasing upwards.
    # Oceananigans z_F is decreasing upwards. z_F[k] is above z_F[k+1].
    # w_field[i,j,k] lives at z_F[k+1].
    # Vertical flux at face k (z_F[k+1]) is w[i,j,k] * WetArea_z[i,j,k]
    # Integration from bottom (k=Nz) upwards (k=Nz to 1).
    # Flux_out_of_cell_k (at face k) = Flux_into_cell_k (at face k-1) + Net_Horiz_Flux_INTO_cell_k

    # Boundary condition: w = 0 at bottom solid boundary (z_F[Nz+1])
    # Flux at bottom face (k=Nz+1) is 0.
    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny
        vertical_volume_flux.data[i+grid.Hx, j+grid.Hy, grid.Nz+grid.Hz] = 0.0 # Flux at z_F[Nz+1] is 0

        # Integrate upwards
        for k in grid.Nz : -1 : 1
            # Net Horizontal Volume Flux INTO cell (i,j,k) (located at C,C,C)
            # Flux IN across west face (i) is u[i,j,k] * wet_face_area_x[i,j,k]
            # Flux OUT across east face (i+1) is u[i+1,j,k] * wet_face_area_x[i+1,j,k]
            # wet_face_area_x is at F,C,C. Indexing aligns with u.
            flux_u_w = cc_params.hFacW[i, j, k] > 1e-10 ? u[i, j, k] * cc_params.wet_face_area_x[i, j, k] : 0.0
            flux_u_e = (i <= grid.Nx && cc_params.hFacW[i+1, j, k] > 1e-10) ? u[i+1, j, k] * cc_params.wet_face_area_x[i+1, j, k] : 0.0
            net_horiz_vol_flux_into_cell = flux_u_w - flux_u_e

            # Vertical volume flux OUT of cell k (across bottom face, k)
            # Flux at face k = Flux at face k+1 + Net Horiz Flux INTO cell k
            # Vertical flux at face k (z_F[k+1]) = vertical_volume_flux[i,j,k]
            # Vertical flux at face k+1 (z_F[k+2]) = vertical_volume_flux[i,j,k+1]

            # Correcting index mapping: w_field[i,j,k] is at z_F[k+1].
            # vertical_volume_flux[i,j,k] should store the flux at z_F[k+1].
            # Integrate upwards from k=Nz to 1.
            # vertical_volume_flux[i,j,k] (Flux at z_F[k+1])
            # Flux_out_of_cell_k = vertical_volume_flux[i,j,k] # flux at z_F[k+1]
            # Flux_into_cell_k = vertical_volume_flux[i,j,k-1] # flux at z_F[k]

            # Let's store the flux at the face index k. w_field[i,j,k] is at face k.
            # Integrate flux starting from k=Nz+1 (bottom, flux=0) upwards to k=1 (top).
            # Vertical flux at face k (w[i,j,k] location)
            # F_k = F_{k+1} + Net_Horiz_Flux_k

            # vertical_volume_flux.data[i+grid.Hx, j+grid.Hy, k+grid.Hz-1] stores flux at face k (z_F[k+1])
            # vertical_volume_flux.data[i+grid.Hx, j+grid.Hy, k+grid.Hz] stores flux at face k+1 (z_F[k+2])

            # Flux at face k (above face k+1) = Flux at face k+1 + Net horizontal flux INTO cell k+1
            # Need to iterate over C cells, not W faces.
            # Loop C cells (i,j,k): 1:Nx, 1:Ny, 1:Nz.
            # Calculate horizontal flux divergence for cell (i,j,k): (u*A)_W - (u*A)_E
            # (u[i,j,k]*wet_face_area_x[i,j,k]) - (u[i+1,j,k]*wet_face_area_x[i+1,j,k])
            # This is volume flux divergence. Volume flux into cell.
            # Sum this divergence from bottom up to get vertical flux.
            # vvf[k] = vvf[k+1] + horiz_div[k] * volume[k] / dz[k]? No.

            # Restart vertical flux calculation:
            # Initialize flux at the bottommost horizontal face (k=Nz+1) to 0.
            # Index in w_field is Nz+1, index in vertical_volume_flux is Nz+1
            w_field.data[i+grid.Hx, j+grid.Hy, grid.Nz+grid.Hz] = 0.0 # w at bottom boundary is 0
            vertical_volume_flux.data[i+grid.Hx, j+grid.Hy, grid.Nz+grid.Hz] = 0.0 # Vertical volume flux at bottom face is 0

            # Integrate upwards from k = grid.Nz to 1 (looping cell indices)
            @inbounds for k_cell in grid.Nz : -1 : 1
                # Horizontal Volume Flux divergence at cell (i,j,k_cell)
                 flux_u_w = cc_params.hFacW[i, j, k_cell] > 1e-10 ? u[i, j, k_cell] * cc_params.wet_face_area_x[i, j, k_cell] : 0.0
                 flux_u_e = (i <= grid.Nx && cc_params.hFacW[i+1, j, k_cell] > 1e-10) ? u[i+1, j, k_cell] * cc_params.wet_face_area_x[i+1, j, k_cell] : 0.0
                 horiz_vol_flux_div_cell = flux_u_e - flux_u_w # Divergence = Out - In

                # Net horizontal volume flux IN = - Divergence
                net_horiz_vol_flux_in_cell = -horiz_vol_flux_div_cell

                # Vertical volume flux at face k_cell (z_F[k_cell+1]) is stored in w_field[i,j,k_cell]
                # Flux_out_of_cell = Flux_into_cell + Net_Horiz_Flux_In_Cell
                # Flux at face k_cell (bottom of cell k_cell) = Flux at face k_cell-1 (top of cell k_cell) + Net Horiz Flux In
                # F_{k_cell+1} = F_{k_cell} + Net_Horiz_Flux_In (for cell k_cell)
                # vertical_volume_flux[i,j,k_cell] = vertical_volume_flux[i,j,k_cell+1] + Net_Horiz_Flux_In_Cell
                # Sum of volume fluxes over boundaries of cell (i,j,k) must be zero.
                # (uA)_E - (uA)_W + (vA)_N - (vA)_S + (wA)_Top - (wA)_Bottom = 0
                # In 2D, v=0. (uA)_E - (uA)_W + (wA)_Top - (wA)_Bottom = 0
                # (wA)_Top = (wA)_Bottom + (uA)_W - (uA)_E
                # (wA) at face k = (wA) at face k+1 + (uA) at face i - (uA) at face i+1
                # Integrate from bottom up (face k=Nz+1 to face k=1)

                # Vertical volume flux at face k_cell (z_F[k_cell+1])
                # Face k_cell is the bottom face of cell (i,j,k_cell).
                # Face k_cell-1 is the top face of cell (i,j,k_cell).
                # Horizontal fluxes into cell (i,j,k_cell) are across faces i and i+1.
                # (u*A)_W = u[i,j,k_cell] * wet_face_area_x[i,j,k_cell] # at Face i
                # (u*A)_E = u[i+1,j,k_cell] * wet_face_area_x[i+1,j,k_cell] # at Face i+1

                 flux_u_w = cc_params.hFacW[i, j, k_cell] > 1e-10 ? u[i, j, k_cell] * cc_params.wet_face_area_x[i, j, k_cell] : 0.0
                 flux_u_e = (i <= grid.Nx && cc_params.hFacW[i+1, j, k_cell] > 1e-10) ? u[i+1, j, k_cell] * cc_face_area_x_at_c_cell_boundaries[i+1, j, k_cell] : 0.0 # Need area at cell boundaries C-grid

                 # This is getting messy with C-grid face areas at C-cell boundaries vs F-points.
                 # Let's go back to the original logic in the provided script which seemed to work.
                 # It computed Net_Horiz_Flux_INTO_cell_k and added it upwards.
                 # This implied (w*A)_k = (w*A)_{k+1} + Net_Horiz_Flux_INTO_cell_k
                 # Vertical volume flux OUT of cell k (across bottom face, k) = Flux_into_cell_k (across top face, k-1) + Net_Horiz_Flux_INTO_cell_k

                 # Flux_out_of_cell_k is the flux across face k (at z_F[k+1]), stored in vertical_volume_flux[i,j,k]
                 # Flux_into_cell_k is the flux across face k-1 (at z_F[k]), stored in vertical_volume_flux[i,j,k-1]

                 # So, vertical_volume_flux[i,j,k_cell] = vertical_volume_flux[i,j,k_cell-1] + Net_Horiz_Flux_INTO_cell_k_cell
                 # Integrating from bottom up:
                 # vertical_volume_flux[i,j,Nz+1] = 0
                 # vertical_volume_flux[i,j,Nz] = vertical_volume_flux[i,j,Nz+1] + Net_Horiz_Flux_INTO_cell_Nz
                 # vertical_volume_flux[i,j,Nz-1] = vertical_volume_flux[i,j,Nz] + Net_Horiz_Flux_INTO_cell_Nz-1
                 # ...
                 # vertical_volume_flux[i,j,k] = vertical_volume_flux[i,j,k+1] + Net_Horiz_Flux_INTO_cell_k

                 # Loop k from Nz down to 1 (C cell index)
                 # Net horizontal flux into cell k: (uA)_W - (uA)_E
                 flux_u_w = cc_params.hFacW[i, j, k_cell] > 1e-10 ? u[i, j, k_cell] * cc_params.wet_face_area_x[i, j, k_cell] : 0.0
                 flux_u_e = (i <= grid.Nx && cc_params.hFacW[i+1, j, k_cell] > 1e-10) ? u[i+1, j, k_cell] * cc_params.wet_face_area_x[i+1, j, k_cell] : 0.0
                 net_horiz_flux_into_cell = flux_u_w - flux_u_e

                 # Vertical volume flux at face k_cell (w_field index k_cell)
                 # This face is the bottom face of cell k_cell.
                 # Flux at face k_cell = Flux at face k_cell + 1 + Net Horiz Flux INTO cell k_cell + 1?
                 # Let's use the provided script's logic:
                 # "Vertical flux across south face (k) ... Flux across north face (k-1)"
                 # This means w[i,j,k] is flux across face k.
                 # Flux at face k = Flux at face k-1 + Net Horiz Flux IN cell k?
                 # This implies integration downwards.

                 # Let's stick to integrating divergence from bottom up.
                 # d(wA)/dz = - d(uA)/dx
                 # Integrate from bottom (z_bottom) to z: w(z)A(z) = w(z_bottom)A(z_bottom) - Integral(d(uA)/dx dz)
                 # w(k)A(k) = w(k+1)A(k+1) - Sum_{k'=k+1 to Nz} (d(uA)/dx)_cell(k') * dz(k')
                 # (d(uA)/dx)_cell(k) = (uA)_E - (uA)_W

                 # Vertical volume flux at face k_cell (w_field index k_cell)
                 # Index k_cell corresponds to z_F[k_cell+1].
                 # w_field[i,j,k_cell] lives at face k_cell.

                 # Flux at face k (w_field index k) = Flux at face k+1 (w_field index k+1) + Horizontal Divergence * Area / dz
                 # Horizontal Divergence (volume flux per unit volume)
                 # div_h_vol = ( (uA)_E - (uA)_W ) / V_cell
                 # V_cell = dx * dy * dz
                 # Horiz Vol Flux Divergence per unit height = (uA)_E - (uA)_W / dz
                 # Horizontal Volume Flux divergence / dy = (uA)_E - (uA)_W / (dy*dz) ?

                # Let's revert to the original continuity diagnostic in the script.
                # It computed a `net_horiz_flux_into_cell` for cell `k_cell` and added it to `vertical_flux_in_k`.
                # It then stored this as `vertical_flux_out_k` at index `k_cell`.
                # This means it was computing: `flux_at_face_k_cell = flux_at_face_k_cell-1 + Net_Horiz_Flux_INTO_cell_k_cell`.
                # This is integrating downwards. Flux at face k-1 = Flux at face k + Net Horiz Flux In cell k.
                # F_{k-1} = F_k + H_k
                # F_0 (top) = F_1 + H_1 = F_2 + H_2 + H_1 = ... = F_{Nz} + Sum(H_{k=1 to Nz})
                # F_{Nz} (bottom face of cell Nz) = F_{Nz+1} (bottom boundary) + Net Horiz Flux IN cell Nz
                # Since F_{Nz+1} = 0, F_{Nz} = Net_Horiz_Flux_INTO_cell_Nz.
                # F_{Nz-1} = F_{Nz} + Net_Horiz_Flux_INTO_cell_Nz-1 = Net_Horiz_Flux_INTO_cell_Nz + Net_Horiz_Flux_INTO_cell_Nz-1
                # So, Flux at face k (w_field index k) = Sum of Net Horiz Flux IN for cells k+1 to Nz.

                # Let's recalculate the horizontal net flux into cell k_cell
                 flux_u_w = cc_params.hFacW[i, j, k_cell] > 1e-10 ? u[i, j, k_cell] * cc_params.wet_face_area_x[i, j, k_cell] : 0.0
                 flux_u_e = (i <= grid.Nx && cc_params.hFacW[i+1, j, k_cell] > 1e-10) ? u[i+1, j, k_cell] * cc_params.wet_face_area_x[i+1, j, k_cell] : 0.0
                 net_horiz_flux_into_cell = flux_u_w - flux_u_e # Net flux INTO cell (i,j,k_cell)

                 # Vertical volume flux at face k_cell (w_field index k_cell)
                 # This flux is Sum of net horizontal fluxes into cells k_cell+1 to Nz.
                 # We need to sum from k_cell+1 down to Nz.
                 # A better way: Integrate divergence from bottom up.
                 # w(z) = w(z_bottom) - Integral(d(uA)/dx dz) / A(z)

                 # Let's try building the horizontal divergence field first.
                 # horiz_vol_flux_div = Field{Center, Center, Center}(grid)
                 # @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
                 #      flux_u_w = cc_params.hFacW[i, j, k] > 1e-10 ? u[i, j, k] * cc_params.wet_face_area_x[i, j, k] : 0.0
                 #      flux_u_e = (i <= grid.Nx && cc_params.hFacW[i+1, j, k] > 1e-10) ? u[i+1, j, k] * cc_params.wet_face_area_x[i+1, j, k] : 0.0
                 #      horiz_vol_flux_div[i, j, k] = (flux_u_e - flux_u_w) # Divergence = Out - In
                 # end
                 # fill_halo_regions!(horiz_vol_flux_div)

                 # Integrate upwards from k=Nz+1 (bottom face)
                 # vertical_volume_flux at face k (w[i,j,k] location)
                 # vvf[i,j,Nz+1] = 0
                 # vvf[i,j,k] = vvf[i,j,k+1] + horiz_vol_flux_div[i,j,k+1] # Flux into cell k+1 at face k+1

                 # The original script's continuity calculation:
                 # Flux_out_of_cell_k (across bottom face k) = Flux_into_cell_k (across top face k-1) + Net_Horiz_Flux_INTO_cell_k
                 # vvf[i,j,k] = vvf[i,j,k-1] + net_horiz_flux_into_cell for cell k

                 # Let's use that logic again, but integrate upwards.
                 # Let vvf[i,j,k] be the vertical volume flux across FACE k (at z_F[k+1])
                 # vvf[i,j,Nz+1] = 0 (Bottom boundary)
                 # For k_cell from Nz down to 1:
                 #   net_horiz_flux_into_cell_k_cell = ...
                 #   vvf[i,j,k_cell] = vvf[i,j,k_cell+1] + net_horiz_flux_into_cell_k_cell

                 # Loop k_cell from Nz down to 1
                 k_face_below = k_cell + 1 # Index of face below cell k_cell (w_field index k_cell+1)
                 k_face_at = k_cell     # Index of face at cell k_cell bottom (w_field index k_cell)

                 # Get flux from face below (already computed in previous iteration)
                 flux_below = vertical_volume_flux.data[i+grid.Hx, j+grid.Hy, k_face_below+grid.Hz-1] # Flux at face k_cell+1

                 # Net horizontal flux into cell k_cell (C index)
                 flux_u_w = cc_params.hFacW[i, j, k_cell] > 1e-10 ? u[i, j, k_cell] * cc_params.wet_face_area_x[i, j, k_cell] : 0.0
                 flux_u_e = (i <= grid.Nx && cc_params.hFacW[i+1, j, k_cell] > 1e-10) ? u[i+1, j, k_cell] * cc_params.wet_face_area_x[i+1, j, k_cell] : 0.0
                 net_horiz_flux_into_cell = flux_u_w - flux_u_e

                 # Vertical volume flux at face k_cell (w_field index k_cell)
                 flux_at_face = flux_below + net_horiz_flux_into_cell
                 vertical_volume_flux.data[i+grid.Hx, j+grid.Hy, k_face_at+grid.Hz-1] = flux_at_face

                 # Convert flux to velocity W at face k_cell
                 if cc_params.wet_face_area_z[i, j, k_face_at] > 1e-15
                      w_field.data[i+grid.Hx, j+grid.Hy, k_face_at+grid.Hz-1] = flux_at_face / cc_params.wet_face_area_z[i, j, k_face_at]
                 else
                      w_field.data[i+grid.Hx, j+grid.Hy, k_face_at+grid.Hz-1] = 0.0
                 end

                 # Ensure W is zero if the face hFacS is zero (explicit mask)
                 if cc_params.hFacS[i, j, k_face_at] < 1e-10
                      w_field.data[i+grid.Hx, j+grid.Hy, k_face_at+grid.Hz-1] = 0.0
                  end

            end # Loop over k_cell (Nz down to 1)
    end # Loop over i,j

    # Handle top boundary w (at k=1, z_F[1]). This w should be zero for a rigid lid.
    # The integration calculated the flux at k=1. If z_F[1] is not the surface, this is an internal face.
    # If z_F[1] is the surface (k=1), then the boundary condition w=0 applies here.
    # The current loop calculates w at k=1 based on flux convergence below.
    # If top BC is w=0, need to enforce it. Let's assume w=0 at z=0 (rigid lid).
    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny
        if grid.zᵃᵃᶠ[1] == 0.0 # If the top face is at z=0
             w_field.data[i+grid.Hx, j+grid.Hy, 1+grid.Hz-1] = 0.0 # Set w at k=1 to 0
         end
    end


    return nothing
end


# --- Model Setup ---

# Set standard advection and closure to 'nothing' to use custom forcings exclusively
advection = nothing_advection()
closure = nothing_closure()

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=alpha, β=beta, T₀=T0, S₀=S0))

# Boundary conditions: Solid wall (zero velocity, zero flux) on sides (x) and bottom (z) by default for Bounded.
# Top (z=0) BC for u: Wind stress applied as a kinematic flux.
wind_stress_flux = wind_stress_tau_x / rho0

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(wind_stress_flux))
# No explicit BCs needed for T and S with custom diffusion forcings if flux = -K*grad is used and grad is handled at boundary.
# For a passive tracer/heat, no flux at solid boundaries is the default zero gradient (∂T/∂n = 0).
# Our custom diffusion handles this by calculating gradient only if both adjacent cells are wet.
# For the top surface, need a BC for T. Let's assume zero flux (insulation).

model = HydrostaticFreeSurfaceModel(
    grid = grid,
    advection = advection, # Standard advection is off
    coriolis = nothing,
    buoyancy = buoyancy,
    closure = closure,     # Standard diffusion/viscosity is off
    boundary_conditions = (u=u_bcs,), # Apply wind stress BC to the u field
    tracers = (:T, :S),
    parameters = (; cut_cell_params), # Pass parameters including hFacs and wet areas/volumes
    forcings = (u=u_forcings, T=T_forcings, S=S_forcings), # Apply all custom forcings
)

# --- Initial Conditions ---

# Linear stratification
T_init(x, y, z) = T0 - 10.0 * (z - grid.zᵃᵃᶠ[1]) / (grid.zᵃᵃᶠ[grid.Nz+1] - grid.zᵃᵃᶠ[1]) # Colder deeper
S_init(x, y, z) = S0 + 0.5 * (z - grid.zᵃᵃᶠ[1]) / (grid.zᵃᵃᶠ[grid.Nz+1] - grid.zᵃᵃᶠ[1]) # Saltier deeper

# Set initial temperature and salinity
set!(model, T=T_init, S=S_init)

# Initialize velocities to zero
set!(model, u=0.0, v=0.0, w=0.0)

# Apply masks to initial conditions based on the *merged* hFacC
# Ensure tracers are zero in dry cells (including those merged away).
@info "Masking initial T and S fields with merged hFacC..."
@inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
    if cut_cell_params.hFacC[i, j, k] < 1e-10 # Use the merged hFacC
        model.tracers.T.data[i+grid.Hx, j+grid.Hy, k+grid.Hz] = 0.0
        model.tracers.S.data[i+grid.Hx, j+grid.Hy, k+grid.Hz] = 0.0
    end
end
# Velocities should also be zero in dry/solid regions.
# u is at Face x. If hFacW=0, u should be 0. This is handled by BCs or model update.
# w is at Face z. If hFacS=0, w should be 0. This is handled by BCs or continuity.
# Let's explicitly zero u where hFacW is zero initially, and w where hFacS is zero.
@info "Masking initial u field with hFacW..."
@inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
     if cut_cell_params.hFacW[i, j, k] < 1e-10
         model.velocities.u.data[i+grid.Hx-1, j+grid.Hy, k+grid.Hz] = 0.0
     end
end
@info "Masking initial w field with hFacS..."
@inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz+1
     if cut_cell_params.hFacS[i, j, k] < 1e-10
         model.velocities.w.data[i+grid.Hx, j+grid.Hy, k+grid.Hz-1] = 0.0
     end
end


# --- Custom Diagnostics ---
# The w diagnostic calculates w *after* the state update.
# This updated w will be used by the advection forcing in the *next* time step.
w_diagnostic = Diagnostic(diagnose_cut_cell_w!, model, field_dependencies = (model.velocities.u,), field_outputs = (w=cut_cell_w,))


# --- Simulation ---
# Time-stepping
simulation = Simulation(model, Δt=dt, stop_time=1000 * dt) # Run for some time

# Add the w diagnostic to the simulation diagnostics
simulation.diagnostics[:cut_cell_w] = w_diagnostic

# Save output
filename = "cut_cell_merged_output_v5.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(
    model,
    (T=model.tracers.T, u=model.velocities.u, S=model.tracers.S, hFacC=model.parameters.cut_cell_params.hFacC,
     wet_cell_volume=model.parameters.cut_cell_params.wet_cell_volume, # Save merged volume
     wet_face_area_x=model.parameters.cut_cell_params.wet_face_area_x, # Save (non-merged) face areas
     wet_face_area_z=model.parameters.cut_cell_params.wet_face_area_z,
     pressure=model.pressure,
    ),
    other_fields = (; w=cut_cell_w), # Output our custom 'w' field
    filename = filename,
    schedule = IterationInterval(50),
    overwrite_existing = true
)

# Progress reporting
progress(sim) = @printf("i: %d, t: %s, max|u|: %.2e, max|w|: %.2e, max|T|: %.2f\n",
                         iteration(sim), prettytime(time(sim)), maximum(abs, sim.model.velocities.u), maximum(abs, sim.diagnostics[:cut_cell_w].field), maximum(abs, sim.model.tracers.T))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

# Run the simulation
println("Starting simulation with cell merging...")
run!(simulation)
println("Simulation complete.")

# --- Analysis and Plotting ---

using JLD2, FieldMetadata
using CairoMakie
# load data
file = jldopen(filename)

iterations = parse.(Int, keys(file["timeseries/t"]))

# --- Create animation ---
anim = @animate for i in iterations

    T = file["timeseries/T/$i"][:, 1, :]
    u = file["timeseries/u/$i"][:, 1, :]
    w = file["timeseries/w/$i"][:, 1, :] # Get the custom diagnosed w
    hFacC = file["timeseries/hFacC/$i"][:, 1, :]
    p = file["timeseries/pressure/$i"][:, 1, :]
    wet_vol = file["timeseries/wet_cell_volume/$i"][:, 1, :] # Get merged volume

    # Extract coordinate arrays
    xC = xnodes(hFacC, Center(), Center(), Center())[:, 1]
    zC = znodes(hFacC, Center(), Center(), Center())[1, 1, :] # z at C
    xF = xnodes(u, Face(), Center(), Center())[:, 1]
    zCu = znodes(u, Face(), Center(), Center())[1, 1, :] # z at u (C)
    zCw = znodes(w, Center(), Center(), Face())[1, 1, :] # z at w (F)

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

    hm_p = heatmap!(ax_p, xC, zC, p, colormap = :viridis)
    Colorbar(fig[2, 2][2, 1], hm_p)

    hm_hFacC = heatmap!(ax_hFacC, xC, zC, hFacC, colormap = :Blues, colorrange=(0, 1))
    Colorbar(fig[3, 1][2, 1], hm_hFacC)

    hm_wetvol = heatmap!(ax_wetvol, xC, zC, wet_vol, colormap = :Plasma)
    Colorbar(fig[3, 2][2, 1], hm_wetvol)

    # Add a title above the plots indicating the time
    Label(fig[0, :], text = "Cut Cell Test with Merging (Flux-Based Forcings): " * string(file["timeseries/t/$i"][end] / day) * " days",
          fontsize = 30)

    # Set axis labels and limits
    for ax in [ax_T, ax_u, ax_w, ax_p, ax_hFacC, ax_wetvol]
        ax.xlabel = "x (m)"
        ax.ylabel = "z (m)"
        ylims!(ax, -total_depth, 0)
        xlims!(ax, 0, total_width)
    end

    fig
end

save("CutCellTest_merged_v5.gif", anim, fps = 10)
close(file)
println("Animation saved in CutCellTest_merged_v5.gif")

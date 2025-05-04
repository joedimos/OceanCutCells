module Forcings

using Oceananigans
using Oceananigans.Grids: ImmersedBoundaryGrid, is_immersed_cell, is_immersed_face,
                          volume, areaᶜᶜᶠ, areaᶠᶜᶜ, # Access geometry from IBG
                          Δx_ᶠᶜᶜ, Δz_ᵃᵃᶠ, Δzᵃᵃᶜ # Access grid spacing functions
using Oceananigans.Operators: ∂xᶠᶜᶜ # Use standard gradient operator
using Oceananigans.Fields: Center, Face, fill_halo_regions
using Oceananigans.Forcings: CustomForcing

# Note: The custom forcing functions now get geometry information from `model.grid`
# and physical parameters from `model.parameters.cut_cell_params`.

# 1. Pressure Gradient Force for u (at Face, Center, Center)
# Adds ∂u∂t = -(1/ρ₀) ∂p/∂x
# Uses Oceananigans' built-in gradient operator ∂xᶠᶜᶜ, which should work on IBG.
# The resulting gradient is then masked by the wet face area/hFacW equivalent.
function add_cut_cell_pressure_gradient_force!(∂u∂t, model)
    grid = model.grid # This is now an ImmersedBoundaryGrid
    p = model.pressure
    cc_params = model.parameters.cut_cell_params # Physical parameters

    # Use Oceananigans' discrete gradient operator ∂xᶠᶜᶜ
    # This operator computes (p[i,j,k] - p[i-1,j,k]) / grid.Δxᶠᶜᶜ[i,j,k]
    # It should handle the IBG correctly, using the appropriate Δx and implicitly zeroing gradient across solid faces.
    p_gradient = ∂xᶠᶜᶜ(1.0 * p) # Compute gradient as a field

    # Add -(1/rho0) * p_gradient to ∂u∂t, masked by wet face/hFacW equivalent
    # Loop over the interior domain of the F,C,C field (∂u∂t)
    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
        # The forcing applies to the u-point (F,C,C).
        # Check if the u-face is wet using is_immersed_face
        if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
             # The gradient operator already handles potential dry points.
             # Add the PGF acceleration term at this wet u-point.
             ∂u∂t[i, j, k] -= (1.0 / cc_params.rho0) * p_gradient[i, j, k]
        else
             # If the u-face is solid, the tendency at this point is effectively zero.
             # Explicitly setting the forcing to zero ensures no contribution.
             ∂u∂t[i, j, k] += 0.0 # Explicitly add zero for clarity
        end
    end
    return nothing
end


# 2. Tracer Advection for T (at Center, Center, Center) - Flux-Based Upwind
# Adds ∂T/∂t = - (1/V_cell) * ∇ . (u T) * V_cell = - ∇ . (u T)
# This uses a manual flux calculation approach. Needs to use IBG geometry accessors.
function add_cut_cell_advection!(∂c∂t, model)
    grid = model.grid # This is now an ImmersedBoundaryGrid
    # Get the tracer field by name from the tendency field's metadata
    c = model.tracers[∂c∂t.name]

    # Need u and w velocities for advection
    u = model.velocities.u
    # w is the diagnosed w field, accessed from model.velocities because we updated model.velocities.w
    w = model.velocities.w

    # cc_params = model.parameters.cut_cell_params # Physical parameters

    # Ensure required fields have halos filled before calculating fluxes
    # (This is typically handled by the time-stepping mechanism, but good to be explicit if needed)
    # fill_halo_regions!(c) # Assuming this is done by Oceananigans before calling forcing
    # fill_halo_regions!(u)
    # fill_halo_regions!(w)


    # Loop over the interior grid points where tracer tendency is computed (C,C,C)
    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
        # Only compute tendency for wet cells using is_immersed_cell
        if !is_immersed_cell(i, j, k, grid, Center(), Center(), Center())

            # --- Horizontal Fluxes (East - West) of tracer c ---
            # Flux across west face (i) = u[i,j,k] * c_upstream * WetArea_face_x[i,j,k]
            # u[i,j,k] is at F[i],C[j],C[k]. Face area is areaᶠᶜᶜ(i,j,k, grid).
            flux_w_c_horiz = 0.0
            # Check if the face is wet using is_immersed_face
            if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
                 # T_upstream: If u[i,j,k] > 0, source is cell (i-1,k). If u[i,j,k] <= 0, source is cell (i,k).
                 # Need to check if source cells are wet using is_immersed_cell.
                 # Be careful with boundary indices (i > 1, i <= grid.Nx)
                 is_wet_im1k = (i > 1) ? !is_immersed_cell(i-1, j, k, grid, Center(), Center(), Center()) : false
                 is_wet_ik = (i <= grid.Nx) ? !is_immersed_cell(i, j, k, grid, Center(), Center(), Center()) : false

                 if is_wet_im1k || is_wet_ik # Flux calculated if face is wet and connects to at least one wet cell
                     u_face = u[i, j, k]
                     if u_face > 0 # Upstream is (i-1, j, k)
                         c_upstream = is_wet_im1k ? c[i-1, j, k] : 0.0 # Use 0 if upstream cell is dry
                     else # Upstream is (i, j, k)
                         c_upstream = is_wet_ik ? c[i, j, k] : 0.0 # Use 0 if upstream cell is dry (should be wet if primary cell is wet)
                     end
                     flux_w_c_horiz = u_face * c_upstream * areaᶠᶜᶜ(i, j, k, grid)
                 end
            end

            # Flux across east face (i+1) = u[i+1,j,k] * c_upstream * WetArea_face_x[i+1,j,k]
            # u[i+1,j,k] is at F[i+1],C[j],C[k]. Face area is areaᶠᶜᶜ(i+1,j,k, grid).
            flux_e_c_horiz = 0.0
            if (i < grid.Nx) && !is_immersed_face(i+1, j, k, grid, Face(), Center(), Center())
                 is_wet_ik = (i <= grid.Nx) ? !is_immersed_cell(i, j, k, grid, Center(), Center(), Center()) : false # Current cell is wet
                 is_wet_ip1k = (i < grid.Nx) ? !is_immersed_cell(i+1, j, k, grid, Center(), Center(), Center()) : false

                 if is_wet_ik || is_wet_ip1k
                    u_face = u[i+1, j, k]
                    if u_face > 0 # Upstream is (i, j, k)
                        c_upstream = is_wet_ik ? c[i, j, k] : 0.0 # Use 0 if upstream cell is dry (should be wet)
                    else # Upstream is (i+1, j, k)
                        c_upstream = is_wet_ip1k ? c[i+1, j, k] : 0.0 # Use 0 if upstream cell is dry
                    end
                    flux_e_c_horiz = u_face * c_upstream * areaᶠᶜᶜ(i+1, j, k, grid)
                 end
            end

            # Horizontal Divergence = (East Flux - West Flux)
            horiz_flux_div_c = flux_e_c_horiz - flux_w_c_horiz

            # --- Vertical Fluxs (North - South) of tracer c ---
            # Flux across south face (k) = w[i,j,k] * c_upstream * WetArea_face_z[i,j,k]
            # w[i,j,k] is at C[i],C[j],F[k]. Face area is areaᶜᶜᶠ(i,j,k, grid).
            # This face is the SOUTHERN face of cell (i,j,k) and NORTHERN face of cell (i,j,k+1).
            flux_s_c_vert = 0.0
             if !is_immersed_face(i, j, k, grid, Center(), Center(), Face())
                  # T_upstream: If w[i,j,k] > 0 (downwards), source is cell (i,k+1). If w[i,j,k] <= 0 (upwards), source is cell (i,k).
                  is_wet_ik = (k <= grid.Nz) ? !is_immersed_cell(i, j, k, grid, Center(), Center(), Center()) : false # Current cell is wet
                  is_wet_ikp1 = (k < grid.Nz) ? !is_immersed_cell(i, j, k+1, grid, Center(), Center(), Center()) : false

                  if is_wet_ik || is_wet_ikp1
                      w_face = w[i, j, k] # w at face k (z_F[k+1])
                      if w_face > 0 # Downwards flux, c from cell (i,k+1)
                          c_upstream = is_wet_ikp1 ? c[i, j, k+1] : 0.0
                      else # Upwards flux, c from cell (i,k)
                          c_upstream = is_wet_ik ? c[i, j, k] : 0.0 # Should be wet
                      end
                      flux_s_c_vert = w_face * c_upstream * areaᶜᶜᶠ(i, j, k, grid)
                  end
             end

            # Flux across north face (k-1) = w[i,j,k-1] * c_upstream * WetArea_face_z[i,j,k-1]
            # w[i,j,k-1] is at C[i],C[j],F[k-1]. Face area is areaᶜᶜᶠ(i,j,k-1, grid).
            # This face is the NORTHERN face of cell (i,j,k) and SOUTHERN face of cell (i,j,k-1).
            flux_n_c_vert = 0.0
            if (k > 1) && !is_immersed_face(i, j, k-1, grid, Center(), Center(), Face())
                  is_wet_ikm1 = (k > 1) ? !is_immersed_cell(i, j, k-1, grid, Center(), Center(), Center()) : false
                  is_wet_ik = (k <= grid.Nz) ? !is_immersed_cell(i, j, k, grid, Center(), Center(), Center()) : false # Current cell is wet

                   if is_wet_ikm1 || is_wet_ik
                       w_face = w[i, j, k-1] # w at face k-1 (z_F[k])
                       if w_face > 0 # Downwards flux, c from cell (i,k)
                           c_upstream = is_wet_ik ? c[i, j, k] : 0.0 # Should be wet
                       else # Upwards flux, c from cell (i,k-1)
                           c_upstream = is_wet_ikm1 ? c[i, j, k-1] : 0.0
                       end
                       flux_n_c_vert = w_face * c_upstream * areaᶜᶜᶠ(i, j, k-1, grid)
                   end
            end

            # Vertical Divergence = (North Flux - South Flux)
            vert_flux_div_c = flux_n_c_vert - flux_s_c_vert

            # Total Divergence
            total_flux_div_c = horiz_flux_div_c + vert_flux_div_c

            # Tendency = - Total Divergence / Wet Cell Volume
            # Wet Cell Volume is volume(i,j,k, grid) for C cell.
            V_cell = volume(i, j, k, grid, Center(), Center(), Center())
            if V_cell > 1e-15 # Avoid division by zero for dry cells
                 ∂c∂t[i, j, k] -= total_flux_div_c / V_cell
            else
                 ∂c∂t[i, j, k] = 0.0 # Tendency is zero in dry cells
            end


        end # if wet cell (!is_immersed_cell)
        # If cell is immersed, ∂c∂t[i,j,k] remains its initialized value (likely 0), which is correct.
    end # Loop
    return nothing
end


# 3. Tracer Diffusion for T and S (at Center, Center, Center) - Flux-Based
# Adds ∂T/∂t = ∇ . (K ∇ T)
function add_cut_cell_diffusion!(∂c∂t, model)
    grid = model.grid # ImmersedBoundaryGrid
    # Get the tracer field by name
    c = model.tracers[∂c∂t.name]
    cc_params = model.parameters.cut_cell_params # Physical parameters
    Kh = cc_params.Kh # Assuming horizontal and vertical diffusivity from params
    Kv = cc_params.Kv

    # fill_halo_regions!(c) # Assuming halos are filled

    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
        # Only compute tendency for wet cells using is_immersed_cell
        if !is_immersed_cell(i, j, k, grid, Center(), Center(), Center())

            # --- Horizontal Fluxes (East - West) of tracer c Gradient ---
            # Flux across west face (i) = -Kh * grad_x(c) * WetArea_face_x
            # Gradient is between c[i-1,j,k] and c[i,j,k]. Distance is Δx_ᶠᶜᶜ(i,j,k, grid).
            # Face area is areaᶠᶜᶜ(i,j,k, grid).
            flux_w_c_horiz = 0.0
            # Check if the face is wet using is_immersed_face
            if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
                 # Gradient requires both adjacent cells to be wet using is_immersed_cell
                 is_wet_im1k = (i > 1) ? !is_immersed_cell(i-1, j, k, grid, Center(), Center(), Center()) : false
                 is_wet_ik = (i <= grid.Nx) ? !is_immersed_cell(i, j, k, grid, Center(), Center(), Center()) : false

                 if is_wet_im1k && is_wet_ik # Gradient calculated only if face connects two wet cells
                      Δx_centers = Δx_ᶠᶜᶜ(i, j, k, grid) # Distance between C points at F location
                      if Δx_centers > 1e-10
                         grad = (c[i, j, k] - c[i-1, j, k]) / Δx_centers
                         flux_w_c_horiz = -Kh * grad * areaᶠᶜᶜ(i, j, k, grid)
                      end
                 end
            end

            # Flux across east face (i+1) = -Kh * grad_x(c) * WetArea_face_x
            # Gradient is between c[i,j,k] and c[i+1,j,k]. Distance is Δx_ᶠᶜᶜ(i+1,j,k, grid).
            # Face area is areaᶠᶜᶜ(i+1,j,k, grid).
            flux_e_c_horiz = 0.0
            if (i < grid.Nx) && !is_immersed_face(i+1, j, k, grid, Face(), Center(), Center())
                 is_wet_ik = (i <= grid.Nx) ? !is_immersed_cell(i, j, k, grid, Center(), Center(), Center()) : false # Current cell is wet
                 is_wet_ip1k = (i < grid.Nx) ? !is_immersed_cell(i+1, j, k, grid, Center(), Center(), Center()) : false

                 if is_wet_ik && is_wet_ip1k
                    Δx_centers = Δx_ᶠᶜᶜ(i+1, j, k, grid)
                    if Δx_centers > 1e-10
                       grad = (c[i+1, j, k] - c[i, j, k]) / Δx_centers
                       flux_e_c_horiz = -Kh * grad * areaᶠᶜᶜ(i+1, j, k, grid)
                    end
                 end
            end

            # Horizontal Flux Divergence = (East Flux - West Flux)
            horiz_flux_div_c = flux_e_c_horiz - flux_w_c_horiz

            # --- Vertical Fluxes (North - South) of tracer c Gradient ---
            # Flux across south face (k) = -Kv * grad_z(c) * WetArea_face_z
            # Gradient is between c[i,j,k] and c[i,j,k+1]. Distance is Δz_ᵃᵃᶠ(i,j,k, grid).
            # Face area is areaᶜᶜᶠ(i,j,k, grid).
            flux_s_c_vert = 0.0
            if !is_immersed_face(i, j, k, grid, Center(), Center(), Face()) # Face must be wet
                 # Gradient requires both adjacent cells to be wet
                 is_wet_ik = (k <= grid.Nz) ? !is_immersed_cell(i, j, k, grid, Center(), Center(), Center()) : false # Current cell is wet
                 is_wet_ikp1 = (k < grid.Nz) ? !is_immersed_cell(i, j, k+1, grid, Center(), Center(), Center()) : false

                 if is_wet_ik && is_wet_ikp1
                      Δz_centers = Δz_ᵃᵃᶠ(i,j,k, grid) # Distance between C points at F location
                      if Δz_centers > 1e-10
                         # c[i,j,k+1] is below c[i,j,k] (z decreases upwards). grad = (c_below - c_above) / dz
                         grad = (c[i, j, k+1] - c[i, j, k]) / Δz_centers
                         # Flux across face k (south face of cell k) is directed downwards if grad is positive (colder above)
                         flux_s_c_vert = -Kv * grad * areaᶜᶜᶠ(i, j, k, grid)
                      end
                 end
            end

            # Flux across north face (k-1) = -Kv * grad_z(c) * WetArea_face_z
            # Gradient is between c[i,j,k-1] and c[i,j,k]. Distance is Δz_ᵃᵃᶠ(i,j,k-1, grid).
            # Face area is areaᶜᶜᶠ(i,j,k-1, grid).
            flux_n_c_vert = 0.0
            if (k > 1) && !is_immersed_face(i, j, k-1, grid, Center(), Center(), Face()) # Face must be wet
                  is_wet_ikm1 = (k > 1) ? !is_immersed_cell(i, j, k-1, grid, Center(), Center(), Center()) : false
                  is_wet_ik = (k <= grid.Nz) ? !is_immersed_cell(i, j, k, grid, Center(), Center(), Center()) : false # Current cell is wet

                  if is_wet_ikm1 && is_wet_ik
                     Δz_centers = Δz_ᵃᵃᶠ(i,j,k-1, grid)
                     if Δz_centers > 1e-10
                         # c[i,j,k-1] is above c[i,j,k]. grad = (c_above - c_below) / dz
                         grad = (c[i, j, k] - c[i, j, k-1]) / Δz_centers # Gradient between k-1 and k
                         # Flux across face k-1 (north face of cell k) is directed upwards if grad is positive (colder below)
                         flux_n_c_vert = -Kv * grad * areaᶜᶜᶠ(i, j, k-1, grid)
                     end
                  end
            end

            # Vertical Flux Divergence = (North Flux - South Flux)
            vert_flux_div_c = flux_n_c_vert - flux_s_c_vert

            # Total Flux Divergence
            total_flux_div_c = horiz_flux_div_c + vert_flux_div_c

            # Tendency = - Total Divergence / Wet Cell Volume
            V_cell = volume(i, j, k, grid, Center(), Center(), Center())
            if V_cell > 1e-15
                 ∂c∂t[i, j, k] -= total_flux_div_c / V_cell
            else
                 ∂c∂t[i, j, k] = 0.0 # Tendency is zero in dry cells
            end

        end # if wet cell (!is_immersed_cell)
    end # Loop
    return nothing
end


# 4. Vertical Advection for u (at Face, Center, Center) - Flux-Based Upwind
# Adds ∂u/∂t = - ∇ . (w u)_vertical
# This calculates ∂/∂z(w u) at F,C,C location.
# Flux is w*u evaluated at (F,C,F) locations. Needs interpolation.
# Approximating w at (C,C,F) and u at (F,C,C) for fluxes at (F,C,F).
function add_cut_cell_vertical_advection_u!(∂u∂t, model)
    grid = model.grid # ImmersedBoundaryGrid
    u = model.velocities.u
    w_field = model.velocities.w # w field updated by diagnostic

    # fill_halo_regions!(u) # Assuming halos are filled
    # fill_halo_regions!(w_field)

    # Loop over the interior domain where u-tendency is calculated (F,C,C)
    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
         # Only compute tendency for wet u-faces using is_immersed_face(F,C,C)
         if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())

            # --- Vertical Fluxes (North - South) of u ---
            # Flux across horizontal face at z_F[k+1] (below u[i,j,k]) - This is face index k.
            # w[i,j,k] is at C[i],C[j],F[k]. WetArea_face_z is areaᶜᶜᶠ(i,j,k, grid).
            # We need flux at F,C,F. Let's approximate w at F by taking w[i,j,k].
            # Flux_below = w[i,j,k] * u_upstream * WetArea_face_z[i,j,k]
            flux_below_u_vert = 0.0 # Flux out the bottom of u[i,j,k] cell
            # Check if the face below u[i,j,k] is wet. This is the face at z_F[k+1], index k.
            if !is_immersed_face(i, j, k, grid, Center(), Center(), Face()) # This face is at C,C,F index k
                 # u_upstream: If w[i,j,k] > 0 (downwards), source is u[i,j,k+1]. If w[i,j,k] <= 0 (upwards), source is u[i,j,k].
                 # Need to check if source u points (F,C,C) are wet.
                 is_wet_u_ik = (k <= grid.Nz) ? !is_immersed_face(i, j, k, grid, Face(), Center(), Center()) : false # u[i,j,k] wet? (Should be true if this loop iteration is reached)
                 is_wet_u_ikp1 = (k < grid.Nz) ? !is_immersed_face(i, j, k+1, grid, Face(), Center(), Center()) : false # u[i,j,k+1] wet?

                 if is_wet_u_ik || is_wet_u_ikp1 # Flux calculated if face is wet and connects to at least one wet u-point
                     # w[i,j,k] is at C,C,F. We need w at F,C,F. Let's use w[i,j,k] as approximation.
                     w_face = w_field[i, j, k] # w velocity at the face (C,C,F index k)
                     # Area at this face is at C,C,F index k. We need F,C,F. Let's use areaᶜᶜᶠ(i,j,k, grid) as approximation.
                     area_face_z = areaᶜᶜᶠ(i, j, k, grid)

                     if w_face > 0 # Downwards flux, upstream is u[i,j,k+1] (at F,C,C index k+1)
                         u_upstream = is_wet_u_ikp1 ? u[i, j, k+1] : 0.0 # If upstream u is dry, use 0
                     else # Upwards flux, upstream is u[i,j,k] (at F,C,C index k)
                         u_upstream = is_wet_u_ik ? u[i, j, k] : 0.0 # Should be wet
                     end
                     flux_below_u_vert = w_face * u_upstream * area_face_z
                 end
            end

            # Flux across horizontal face at z_F[k] (above u[i,j,k]) - This is face index k-1.
            # w[i,j,k-1] is at C[i],C[j],F[k-1]. WetArea_face_z is areaᶜᶜᶠ(i,j,k-1, grid).
            # Flux_above = w[i,j,k-1] * u_upstream * WetArea_face_z[i,j,k-1]
            flux_above_u_vert = 0.0 # Flux into the top of u[i,j,k] cell
            # Check if the face above u[i,j,k] is wet. This is the face at z_F[k], index k-1.
            if (k > 1) && !is_immersed_face(i, j, k-1, grid, Center(), Center(), Face()) # Face is at C,C,F index k-1

                  is_wet_u_ikm1 = (k > 1) ? !is_immersed_face(i, j, k-1, grid, Face(), Center(), Center()) : false # u[i,j,k-1] wet?
                  is_wet_u_ik = (k <= grid.Nz) ? !is_immersed_face(i, j, k, grid, Face(), Center(), Center()) : false # u[i,j,k] wet? (Should be true)

                   if is_wet_u_ikm1 || is_wet_u_ik # Flux calculated if face is wet and connects to at least one wet u-point
                      w_face = w_field[i, j, k-1] # w velocity at the face (C,C,F index k-1)
                      area_face_z = areaᶜᶜᶠ(i, j, k-1, grid) # Area at this face (C,C,F index k-1)

                       if w_face > 0 # Downwards flux, upstream is u[i,j,k]
                           u_upstream = is_wet_u_ik ? u[i, j, k] : 0.0 # Should be wet
                       else # Upwards flux, upstream is u[i,j,k-1]
                           u_upstream = is_wet_u_ikm1 ? u[i, j, k-1] : 0.0
                       end
                       flux_above_u_vert = w_face * u_upstream * area_face_z
                   end
            end

            # Vertical Divergence = (Flux_into_top - Flux_out_bottom) / Volume_u_cell
            # Volume of u cell is volume(i,j,k, grid) at F,C,C location.
            V_u_cell = volume(i, j, k, grid, Face(), Center(), Center()) # Use IBG volume function
            if V_u_cell > 1e-15
                 # Divergence of vertical flux = (Flux_above - Flux_below) / V_u_cell * dz? No, Volume.
                 # div = (Flux_above - Flux_below) / V_u_cell * Area_face_z? No.
                 # div = (Flux_above - Flux_below) / Δz_u_cell ? Yes, for acceleration units.
                 # Δz_u_cell = grid.Δzᵃᵃᶜ[i,j,k] # Use grid accessor for u point

                 # The divergence is (flux_above - flux_below) / V_u_cell * Area_face_z? No.
                 # The divergence is d(w u)/dz. Discretely this is [(wu)_above - (wu)_below] / Δz_u_cell
                 # We computed flux = wu * Area_face. So we need (Flux_above/Area_above - Flux_below/Area_below) / Δz_u_cell?
                 # No, the tendency is - ∇ . (uv). For vertical part: -∂/∂z(wu).
                 # Discrete form at F,C,C (i,j,k) is - [ (wu)_face_k-1 - (wu)_face_k ] / Δz_u_cell
                 # where (wu)_face_k-1 = Flux_above / Area_face_above
                 #       (wu)_face_k = Flux_below / Area_face_below

                 # This is getting complex with staggered grid and cut cells.
                 # Let's approximate the vertical divergence simply as:
                 # (Flux_above - Flux_below) / V_u_cell (Volume divergence -> m^3/s / m^3 = 1/s)
                 # We need acceleration (m/s^2). The forcings are defined to add ACCELERATION.
                 # ∇ . (uv) has units of velocity * gradient(velocity) ~ (m/s) * ( (m/s)/m ) = m/s^2.
                 # So the formula should be -∂/∂z(wu) which is approx - [(wu)_above - (wu)_below] / Δz_u_cell.
                 # where (wu)_above and (wu)_below are velocities at the faces.
                 # (wu)_face_k-1 is the vertical velocity at face k-1 multiplied by u interpolated at Face k-1.
                 # Let's use the fluxes calculated:
                 # vertical_advection_accel = - (Flux_above_u_vert - Flux_below_u_vert) / V_u_cell * Area_face_z ? No.

                 # Let's use the simpler finite volume form: Sum of fluxes / Volume
                 # ∂u/∂t = - (Flux_out - Flux_in) / Volume = - (Flux_below - Flux_above) / V_u_cell
                 vertical_advection_accel = - (flux_below_u_vert - flux_above_u_vert) / V_u_cell

                 ∂u∂t[i, j, k] += vertical_advection_accel

            else # If the u-face is immersed
                 ∂u∂t[i, j, k] = 0.0 # Tendency is zero
             end
        end # Loop
        return nothing
    end


# 5. Vertical Diffusion for u (at Face, Center, Center) - Flux-Based Vertical Viscosity
# Adds ∂u/∂t = ∂/∂z (Av ∂u/∂z)
# Uses a centered stencil approximation for the Laplacian, masked by wet face connectivity.
function add_cut_cell_vertical_diffusion_u!(∂u∂t, model)
    grid = model.grid # ImmersedBoundaryGrid
    u = model.velocities.u
    cc_params = model.parameters.cut_cell_params
    Av = cc_params.Av

    # fill_halo_regions!(u) # Assuming halos filled

    # Loop over the interior domain where u-tendency is calculated (F,C,C)
    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
         # Only compute tendency for wet u-faces using is_immersed_face(F,C,C)
         if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())

            # Centered difference approximation for ∂/∂z (Av ∂u/∂z) at (i,j,k)
            # The term is Av * ∂^2 u / ∂z^2, approx Av * [ (u_k+1 - u_k)/dz_k+1/2 - (u_k - u_k-1)/dz_k-1/2 ] / dz_k
            # u_k = u[i,j,k]
            # dz_k+1/2 is distance between u[i,j,k] and u[i,j,k+1], which is Δz at face k (C,C,F index k) = Δz_ᵃᵃᶠ(i,j,k, grid)
            # dz_k-1/2 is distance between u[i,j,k-1] and u[i,j,k], which is Δz at face k-1 (C,C,F index k-1) = Δz_ᵃᵃᶠ(i,j,k-1, grid)
            # dz_k is the height of the u cell (F,C,C index k) = Δzᵃᵃᶜ(i,j,k, grid)

            # Need to check connectivity using is_immersed_face for u points
            is_wet_u_ik = (k <= grid.Nz) ? !is_immersed_face(i, j, k, grid, Face(), Center(), Center()) : false # u[i,j,k] wet? (Should be true)
            is_wet_u_ikm1 = (k > 1) ? !is_immersed_face(i, j, k-1, grid, Face(), Center(), Center()) : false # u[i,j,k-1] wet?
            is_wet_u_ikp1 = (k < grid.Nz) ? !is_immersed_face(i, j, k+1, grid, Face(), Center(), Center()) : false # u[i,j,k+1] wet?

            # Vertical gradient below u[i,j,k] (between k and k+1)
            grad_below = 0.0
            # Need both u points wet AND the connecting face below (at C,C,F index k) to be wet for diffusion
            if is_wet_u_ik && is_wet_u_ikp1 && !is_immersed_face(i, j, k, grid, Center(), Center(), Face())
                Δz_between = Δz_ᵃᵃᶠ(i,j,k, grid) # Distance between u points at F,C,C k and k+1
                 if Δz_between > 1e-10
                    grad_below = (u[i, j, k+1] - u[i, j, k]) / Δz_between
                 end
            end

            # Vertical gradient above u[i,j,k] (between k-1 and k)
            grad_above = 0.0
            # Need both u points wet AND the connecting face above (at C,C,F index k-1) to be wet
             if is_wet_u_ikm1 && is_wet_u_ik && (k > 1) && !is_immersed_face(i, j, k-1, grid, Center(), Center(), Face())
                 Δz_between = Δz_ᵃᵃᶠ(i,j,k-1, grid) # Distance between u points at F,C,C k-1 and k
                 if Δz_between > 1e-10
                     grad_above = (u[i, j, k] - u[i, j, k-1]) / Δz_between
                 end
             end

            # Vertical Laplacian: ∂/∂z(Av ∂u/∂z) approx Av * (grad_above - grad_below) / Δz_u_cell
            Δz_u_cell = Δzᵃᵃᶜ(i, j, k, grid) # Height of the u-cell control volume (F,C,C location)
            if Δz_u_cell > 1e-10
                 laplacian_z = (grad_above - grad_below) / Δz_u_cell
                 ∂u∂t[i, j, k] += Av * laplacian_z
            end

         else # If the u-face is immersed
             ∂u∂t[i, j, k] = 0.0 # Tendency is zero
         end
    end # Loop
    return nothing
end


# 6. Bottom Drag (at Face, Center, Center)
# Adds ∂u/∂t = - Cd |u| u / Δz_bottom_cell
function add_cut_cell_bottom_drag!(∂u∂t, model)
    grid = model.grid # ImmersedBoundaryGrid
    u = model.velocities.u
    cc_params = model.parameters.cut_cell_params # Physical parameters
    Cd = cc_params.bottom_drag_coeff

    # fill_halo_regions!(u) # Assuming halos filled

    # Apply drag to the bottom-most wet u-point in each column
    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny
        k_bottom = -1
        # Find the deepest wet u-face point [i,j,k]. Iterate from bottom up.
        # Use is_immersed_face(F,C,C) to find wet u-points.
        for k in grid.Nz : -1 : 1
           if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
               k_bottom = k
               break # Found the bottom-most wet u point in this column (i,j)
           end
        end

        if k_bottom != -1
            u_bottom = u[i, j, k_bottom]
            # Drag force per unit mass = - Cd * |u_bottom| * u_bottom / H_eff
            # H_eff is the vertical extent over which the drag is distributed. Use the height of the u-cell: Δzᵃᵃᶜ[i,j,k_bottom].
            H_eff = Δzᵃᵃᶜ(i, j, k_bottom, grid) # Use IBG grid accessor

            if H_eff > 1e-10
                drag_accel = - Cd * abs(u_bottom) * u_bottom / H_eff
                ∂u∂t[i, j, k_bottom] += drag_accel
            end
        end
    end
    return nothing
end


# Bundle forcings using the generic advection and diffusion functions
# Note: field_dependencies correctly listed for the generic functions
u_forcings = (
    pressure_gradient = CustomForcing(add_cut_cell_pressure_gradient_force!, field_dependencies=(:p,)), # PGF depends on pressure
    vertical_advection = CustomForcing(add_cut_cell_vertical_advection_u!, field_dependencies=(:u, :w)), # V-Adv depends on u, w
    vertical_diffusion = CustomForcing(add_cut_cell_vertical_diffusion_u!, field_dependencies=(:u,)),   # V-Diff depends on u
    bottom_drag = CustomForcing(add_cut_cell_bottom_drag!, field_dependencies=(:u,)) # Bottom Drag depends on u
)

T_forcings = (
     advection = CustomForcing(add_cut_cell_advection!, field_dependencies=(:T, :u, :w)), # Advection depends on T, u, w
     diffusion = CustomForcing(add_cut_cell_diffusion!, field_dependencies=(:T,)) # Diffusion depends on T
)

S_forcings = (
     advection = CustomForcing(add_cut_cell_advection!, field_dependencies=(:S, :u, :w)), # Advection depends on S, u, w
     diffusion = CustomForcing(add_cut_cell_diffusion!, field_dependencies=(:S,)) # Diffusion depends on S
)


end # module Forcings

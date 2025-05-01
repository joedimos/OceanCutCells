    module Forcings

    using Oceananigans
    using Oceananigans.Grids: Δxᶠᶜᶜ, Δzᵃᵃᶠ # Added grid spacing accessors
    using Oceananigans.Operators: ∂xᶠᶜᶜ
    using Oceananigans.Fields: Center, Face

    using ..Parameters: CutCellParameters # Access the struct
    # Need access to velocities and tracers from the model
    # CustomForcing functions receive (∂c∂t, model)

    # 1. Pressure Gradient Force for u (at Face, Center, Center)
    # Adds ∂u∂t = -(1/ρ₀) ∂p/∂x
    function add_cut_cell_pressure_gradient_force!(∂u∂t, model)
        grid = model.grid
        p = model.pressure
        cc_params = model.parameters.cut_cell_params
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        # Use Oceananigans' discrete gradient operator ∂xᶠᶜᶜ
        # This operator is defined over the internal Face points (1:Nx+1 in x)
        # It computes (p[i,j,k] - p[i-1,j,k]) / grid.Δxᶠᶜᶜ[i,j,k]
        # Note: ∂xᶠᶜᶜ returns a field at F,C,C location. It handles halo regions internally if p has correct halos.
        p_gradient = ∂xᶠᶜᶜ(1.0 * p) # Compute gradient field

        # Add -(1/rho0) * p_gradient to ∂u∂t, masked by hFacW
        @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
            # hFacW gates the flux across the face. If the face is solid (hFacW=0), the PGF should be balanced by the wall.
            # If the face is wet, PGF contributes to acceleration.
            # The PGF is applied at the u-point (F,C,C). The hFacW indicates if this point is active.
            if cc_params.hFacW[i, j, k] > 1e-10
                 ∂u∂t.data[i+Hx-1, j+Hy, k+Hz] -= (1.0 / cc_params.rho0) * p_gradient[i, j, k]
            else
                 # If hFacW is zero, the tendency at this u-point is effectively zero.
                 # Explicitly setting the forcing to zero here ensures no contribution.
                 ∂u∂t.data[i+Hx-1, j+Hy, k+Hz] += 0.0 # Explicitly add zero for clarity
            end
        end
        return nothing
    end


    # 2. Tracer Advection for T (at Center, Center, Center) - Flux-Based Upwind
    # Adds ∂T/∂t = - (1/V_cell) * ∇ . (u T) * V_cell = - ∇ . (u T)
    function add_cut_cell_advection_T!(∂T∂t, model)
        grid = model.grid
        T = model.tracers.T
        u = model.velocities.u # Need u for horizontal advection
        w = model.velocities.w # Need w for vertical advection (this is the DIAGNOSED w field)

        cc_params = model.parameters.cut_cell_params
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        # Ensure required fields have halos filled
        fill_halo_regions!(T)
        fill_halo_regions!(u)
        fill_halo_regions!(w)

        # Loop over the interior grid points where T tendency is computed (C,C,C)
        @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
            # Only compute tendency for wet cells (using the merged hFacC mask)
            if cc_params.hFacC[i, j, k] > 1e-10

                # --- Horizontal Fluxes (East - West) of T ---
                # Flux across west face (i) = u[i,j,k] * T_upstream * WetArea_face_x[i,j,k]
                flux_w_T_horiz = 0.0
                if cc_params.hFacW[i, j, k] > 1e-10
                     # T_upstream: If u[i,j,k] > 0, source is cell (i-1,k). If u[i,j,k] <= 0, source is cell (i,k).
                     is_wet_im1k = (i > 1 && cc_params.hFacC[i-1, j, k] > 1e-10)
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet

                     if is_wet_im1k || is_wet_ik # Flux calculated if face is wet and connects to at least one wet cell
                         u_face = u[i, j, k]
                         if u_face > 0 # Upstream is (i-1, j, k)
                             T_upstream = is_wet_im1k ? T[i-1, j, k] : 0.0 # Use 0 if upstream cell is dry
                         else # Upstream is (i, j, k)
                             T_upstream = is_wet_ik ? T[i, j, k] : 0.0 # Use 0 if upstream cell is dry (should be wet here if cc_params.hFacC[i,j,k] > 1e-10)
                         end
                         flux_w_T_horiz = u_face * T_upstream * cc_params.wet_face_area_x[i, j, k]
                     end
                end

                # Flux across east face (i+1) = u[i+1,j,k] * T_upstream * WetArea_face_x[i+1,j,k]
                flux_e_T_horiz = 0.0
                if cc_params.hFacW[i+1, j, k] > 1e-10
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet
                     is_wet_ip1k = (i < grid.Nx && cc_params.hFacC[i+1, j, k] > 1e-10)

                     if is_wet_ik || is_wet_ip1k
                        u_face = u[i+1, j, k]
                        if u_face > 0 # Upstream is (i, j, k)
                            T_upstream = is_wet_ik ? T[i, j, k] : 0.0 # Use 0 if upstream cell is dry (should be wet)
                        else # Upstream is (i+1, j, k)
                            T_upstream = is_wet_ip1k ? T[i+1, j, k] : 0.0 # Use 0 if upstream cell is dry
                        end
                        flux_e_T_horiz = u_face * T_upstream * cc_params.wet_face_area_x[i+1, j, k]
                     end
                end

                # Horizontal Divergence = (East Flux - West Flux)
                horiz_flux_div_T = flux_e_T_horiz - flux_w_T_horiz

                # --- Vertical Fluxs (North - South) of T ---
                # Flux across south face (k) = w[i,j,k] * T_upstream * WetArea_face_z[i,j,k]
                # Note: w[i,j,k] lives at z_F[k+1] (Face index k+1).
                # wet_face_area_z[i,j,k] is also at Face z index k.
                # There's a mismatch in indexing w and wet_face_area_z in the original code.
                # w[i,j,k] lives at z_F[k+1]. The face between cell k and k+1 is face k+1.
                # hFacS[i,j,k] is at Face z index k, which is z_F[k+1].
                # Let's assume w[i,j,k] corresponds to the vertical velocity at face k.
                # Flux across face k = w[i,j,k] * T_upstream * WetArea_face_z[i,j,k]
                # This face k is the SOUTHERN face of cell (i,j,k) and NORTHERN face of cell (i,j,k+1).

                # Flux across south face (k) of cell (i,j,k) (This is at z_F[k+1]) = w[i,j,k] * T_upstream * WetArea_face_z[i,j,k]
                flux_s_T_vert = 0.0
                 if cc_params.hFacS[i, j, k] > 1e-10
                      # T_upstream: If w[i,j,k] > 0 (downwards), source is cell (i,k+1). If w[i,j,k] <= 0 (upwards), source is cell (i,k).
                      # w[i,j,k] is at z_F[k+1] (Face index k).
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet
                      is_wet_ikp1 = (k < grid.Nz && cc_params.hFacC[i, j, k+1] > 1e-10)

                      if is_wet_ik || is_wet_ikp1
                          w_face = w[i, j, k] # w at face k (z_F[k+1])
                          if w_face > 0 # Downwards flux, T from cell (i,k+1)
                              T_upstream = is_wet_ikp1 ? T[i, j, k+1] : 0.0
                          else # Upwards flux, T from cell (i,k)
                              T_upstream = is_wet_ik ? T[i, j, k] : 0.0 # Should be wet
                          end
                          flux_s_T_vert = w_face * T_upstream * cc_params.wet_face_area_z[i, j, k]
                      end
                 end

                # Flux across north face (k-1) of cell (i,j,k) (This is at z_F[k]) = w[i,j,k-1] * T_upstream * WetArea_face_z[i,j,k-1]
                flux_n_T_vert = 0.0
                if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10
                      is_wet_ikm1 = (k > 1 && cc_params.hFacC[i, j, k-1] > 1e-10)
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet

                       if is_wet_ikm1 || is_wet_ik
                           w_face = w[i, j, k-1] # w at face k-1 (z_F[k])
                           if w_face > 0 # Downwards flux, T from cell (i,k)
                               T_upstream = is_wet_ik ? T[i, j, k] : 0.0 # Should be wet
                           else # Upwards flux, T from cell (i,k-1)
                               T_upstream = is_wet_ikm1 ? T[i, j, k-1] : 0.0
                           end
                           flux_n_T_vert = w_face * T_upstream * cc_params.wet_face_area_z[i, j, k-1]
                       end
                end

                # Vertical Divergence = (Flux_into_top - Flux_out_bottom) = (North Flux - South Flux)
                vert_flux_div_T = flux_n_T_vert - flux_s_T_vert

                # Total Divergence
                total_flux_div_T = horiz_flux_div_T + vert_flux_div_T

                # Tendency = - Total Divergence / Wet Cell Volume (Use merged volume)
                ∂T∂t.data[i+Hx, j+Hy, k+Hz] -= total_flux_div_T * cc_params.recip_wet_cell_volume[i, j, k]

            end # if wet cell (using merged hFacC)
        end # Loop
        return nothing
    end


    # 3. Tracer Diffusion for T (at Center, Center, Center) - Flux-Based
    # Adds ∂T/∂t = ∇ . (K ∇ T)
    function add_cut_cell_diffusion_T!(∂T∂t, model)
        grid = model.grid
        T = model.tracers.T
        cc_params = model.parameters.cut_cell_params
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        fill_halo_regions!(T) # Ensure halos are up-to-date for gradient calculation

        @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
            # Only compute tendency for wet cells (using merged hFacC)
            if cc_params.hFacC[i, j, k] > 1e-10

                # --- Horizontal Fluxes (East - West) of T Gradient ---
                # Flux across west face (i) = -Kh * grad_x(T) * WetArea_face_x[i,j,k]
                flux_w_T_horiz = 0.0
                if cc_params.hFacW[i, j, k] > 1e-10 # Face must be wet
                     # Gradient requires both adjacent C cells to be wet (using merged hFacC)
                     is_wet_im1k = (i > 1 && cc_params.hFacC[i-1, j, k] > 1e-10)
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet

                     if is_wet_im1k && is_wet_ik # Gradient calculated only if face connects two wet cells
                          # Δx_centers = grid.Δxᶠᶜᶜ[i, j, k] # Distance between T at (i-1,k) and (i,k)
                          # Using grid accessor specific to C-cell distance at F location
                          Δx_centers = grid.Δx_ᶠᶜᶜ(i,j,k)
                          if Δx_centers > 1e-10
                             grad = (T[i, j, k] - T[i-1, j, k]) / Δx_centers
                             flux_w_T_horiz = -cc_params.Kh * grad * cc_params.wet_face_area_x[i, j, k]
                          end
                     end
                end

                # Flux across east face (i+1) = -Kh * grad_x(T) * WetArea_face_x[i+1,j,k]
                flux_e_T_horiz = 0.0
                if cc_params.hFacW[i+1, j, k] > 1e-10
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet
                     is_wet_ip1k = (i < grid.Nx && cc_params.hFacC[i+1, j, k] > 1e-10)

                     if is_wet_ik && is_wet_ip1k
                        Δx_centers = grid.Δx_ᶠᶜᶜ(i+1,j,k)
                        if Δx_centers > 1e-10
                           grad = (T[i+1, j, k] - T[i, j, k]) / Δx_centers
                           flux_e_T_horiz = -cc_params.Kh * grad * cc_params.wet_face_area_x[i+1, j, k]
                        end
                     end
                end

                # Horizontal Flux Divergence = (East Flux - West Flux)
                horiz_flux_div_T = flux_e_T_horiz - flux_w_T_horiz

                # --- Vertical Fluxes (North - South) of T Gradient ---
                # Flux across south face (k) = -Kv * grad_z(T) * WetArea_face_z[i,j,k]
                # Flux across face k (at z_F[k+1]). Gradient is between T[i,j,k] and T[i,j,k+1].
                flux_s_T_vert = 0.0
                if cc_params.hFacS[i, j, k] > 1e-10 # Face must be wet
                     # Gradient requires both adjacent C cells to be wet (using merged hFacC)
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet
                     is_wet_ikp1 = (k < grid.Nz && cc_params.hFacC[i, j, k+1] > 1e-10)

                     if is_wet_ik && is_wet_ikp1
                          # Δz_centers = grid.Δzᵃᵃᶠ[k, j, i] # Original had this index order, check grid accessors
                          # Using grid accessor specific to C-cell distance at F location
                          Δz_centers = grid.Δz_ᵃᵃᶠ(i,j,k)
                          if Δz_centers > 1e-10
                             # T[i,j,k+1] is below T[i,j,k] (z decreases upwards). grad = (T_below - T_above) / dz
                             grad = (T[i, j, k+1] - T[i, j, k]) / Δz_centers
                             # Flux across face k (south face of cell k) is directed downwards if grad is positive (colder above)
                             flux_s_T_vert = -cc_params.Kv * grad * cc_params.wet_face_area_z[i, j, k]
                          end
                     end
                end

                # Flux across north face (k-1) = -Kv * grad_z(T) * WetArea_face_z[i,j,k-1]
                # Flux across face k-1 (at z_F[k]). Gradient is between T[i,j,k-1] and T[i,j,k].
                flux_n_T_vert = 0.0
                if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10 # Face must be wet
                      is_wet_ikm1 = (k > 1 && cc_params.hFacC[i, j, k-1] > 1e-10)
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet

                      if is_wet_ikm1 && is_wet_ik
                         Δz_centers = grid.Δz_ᵃᵃᶠ(i,j,k-1)
                         if Δz_centers > 1e-10
                             # T[i,j,k-1] is above T[i,j,k]. grad = (T_above - T_below) / dz
                             grad = (T[i, j, k] - T[i, j, k-1]) / Δz_centers # Gradient between k-1 and k
                             # Flux across face k-1 (north face of cell k) is directed upwards if grad is positive (colder below)
                             flux_n_T_vert = -cc_params.Kv * grad * cc_params.wet_face_area_z[i, j, k-1]
                         end
                      end
                end

                # Vertical Flux Divergence = (North Flux - South Flux)
                vert_flux_div_T = flux_n_T_vert - flux_s_T_vert

                # Total Flux Divergence
                total_flux_div_T = horiz_flux_div_T + vert_flux_div_T

                # Tendency = - Total Divergence / Wet Cell Volume (Use merged volume)
                ∂T∂t.data[i+Hx, j+Hy, k+Hz] -= total_flux_div_T * cc_params.recip_wet_cell_volume[i, j, k]

            end # if wet cell (using merged hFacC)
        end # Loop
        return nothing
    end

    # Custom forcing for S diffusion (same as T diffusion)
    function add_cut_cell_diffusion_S!(∂S∂t, model)
        grid = model.grid
        S = model.tracers.S
        cc_params = model.parameters.cut_cell_params
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        fill_halo_regions!(S) # Ensure halos are up-to-date

        @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
            if cc_params.hFacC[i, j, k] > 1e-10

                # Horizontal Diffusion (same logic as T)
                flux_w_S_horiz = 0.0
                if cc_params.hFacW[i, j, k] > 1e-10
                     is_wet_im1k = (i > 1 && cc_params.hFacC[i-1, j, k] > 1e-10)
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                     if is_wet_im1k && is_wet_ik
                          Δx_centers = grid.Δx_ᶠᶜᶜ(i,j,k)
                          if Δx_centers > 1e-10
                             grad = (S[i, j, k] - S[i-1, j, k]) / Δx_centers
                             flux_w_S_horiz = -cc_params.Kh * grad * cc_params.wet_face_area_x[i, j, k]
                          end
                     end
                end

                flux_e_S_horiz = 0.0
                if cc_params.hFacW[i+1, j, k] > 1e-10
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                     is_wet_ip1k = (i < grid.Nx && cc_params.hFacC[i+1, j, k] > 1e-10)
                     if is_wet_ik && is_wet_ip1k
                        Δx_centers = grid.Δx_ᶠᶜᶜ(i+1,j,k)
                        if Δx_centers > 1e-10
                           grad = (S[i+1, j, k] - S[i, j, k]) / Δx_centers
                           flux_e_S_horiz = -cc_params.Kh * grad * cc_params.wet_face_area_x[i+1, j, k]
                        end
                     end
                end
                horiz_flux_div_S = flux_e_S_horiz - flux_w_S_horiz

                # Vertical Diffusion (same logic as T)
                flux_s_S_vert = 0.0
                 if cc_params.hFacS[i, j, k] > 1e-10
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                      is_wet_ikp1 = (k < grid.Nz && cc_params.hFacC[i, j, k+1] > 1e-10)
                      if is_wet_ik && is_wet_ikp1
                           Δz_centers = grid.Δz_ᵃᵃᶠ(i,j,k)
                           if Δz_centers > 1e-10
                              grad = (S[i, j, k+1] - S[i, j, k]) / Δz_centers
                              flux_s_S_vert = -cc_params.Kv * grad * cc_params.wet_face_area_z[i, j, k]
                           end
                      end
                 end

                flux_n_S_vert = 0.0
                if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10
                      is_wet_ikm1 = (k > 1 && cc_params.hFacC[i, j, k-1] > 1e-10)
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                      if is_wet_ikm1 && is_wet_ik
                         Δz_centers = grid.Δz_ᵃᵃᶠ(i,j,k-1)
                         if Δz_centers > 1e-10
                             grad = (S[i, j, k] - S[i, j, k-1]) / Δz_centers
                             flux_n_S_vert = -cc_params.Kv * grad * cc_params.wet_face_area_z[i, j, k-1]
                         end
                      end
                end
                vert_flux_div_S = flux_n_S_vert - flux_s_S_vert

                # Total Flux Divergence
                total_flux_div_S = horiz_flux_div_S + vert_flux_div_S

                # Tendency = - Total Divergence / Wet Cell Volume (Use merged volume)
                ∂S∂t.data[i+Hx, j+Hy, k+Hz] -= total_flux_div_S * cc_params.recip_wet_cell_volume[i, j, k]

            end
        end
        return nothing
    end


    # 4. Vertical Advection for u (at Face, Center, Center) - Flux-Based Upwind
    # Adds ∂u/∂t = - ∇ . (w u)_vertical
    function add_cut_cell_vertical_advection_u!(∂u∂t, model)
        grid = model.grid
        u = model.velocities.u
        w = model.velocities.w # w field updated by diagnostic
        cc_params = model.parameters.cut_cell_params
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        fill_halo_regions!(u) # Need u halos for upwinding
        fill_halo_regions!(w) # Need w halos

        # Loop over the interior domain where u-tendency is calculated (F,C,C)
        @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
             # Only compute tendency for wet u-faces (using original hFacW)
             if cc_params.hFacW[i, j, k] > 1e-10

                # --- Vertical Fluxes (North - South) of u ---
                # Flux across horizontal face k (at z_F[k+1]) - where w[i,j,k] lives and wet_face_area_z[i,j,k] is defined.
                # This face is the SOUTHERN face of the u-point control volume at (i,j,k).
                # Flux_s = w[i,j,k] * u_upstream * WetArea_face_z[i,j,k]
                flux_s_u_vert = 0.0 # Flux out the bottom of u[i,j,k] cell
                if cc_params.hFacS[i, j, k] > 1e-10 # Face must be wet

                     # Check connectivity of adjacent u-points for upwinding u using original hFacW
                     is_wet_u_ik = cc_params.hFacW[i, j, k] > 1e-10 # u[i,j,k] wet? (Should be true)
                     is_wet_u_ikp1 = (k < grid.Nz && cc_params.hFacW[i, j, k+1] > 1e-10) # u[i,j,k+1] wet?

                     if is_wet_u_ik || is_wet_u_ikp1 # Flux calculated if face is wet and connects to at least one wet u-point
                         w_face = w[i, j, k] # w velocity at face k (z_F[k+1])
                         if w_face > 0 # Downwards flux, upstream is u[i,j,k+1]
                             u_upstream = is_wet_u_ikp1 ? u[i, j, k+1] : 0.0 # If upstream u is dry, use 0
                         else # Upwards flux, upstream is u[i,j,k]
                             u_upstream = is_wet_u_ik ? u[i, j, k] : 0.0 # Should be wet
                         end
                         flux_s_u_vert = w_face * u_upstream * cc_params.wet_face_area_z[i, j, k]
                     end
                end

                # Flux across horizontal face k-1 (at z_F[k]) - where w[i,j,k-1] lives and wet_face_area_z[i,j,k-1] is defined.
                # This face is the NORTHERN face of the u-point control volume at (i,j,k).
                # Flux_n = w[i,j,k-1] * u_upstream * WetArea_face_z[i,j,k-1]
                flux_n_u_vert = 0.0 # Flux into the top of u[i,j,k] cell
                if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10 # Face must be wet

                      is_wet_u_ikm1 = (k > 1 && cc_params.hFacW[i, j, k-1] > 1e-10) # u[i,j,k-1] wet?
                      is_wet_u_ik = cc_params.hFacW[i, j, k] > 1e-10 # u[i,j,k] wet? (Should be true)

                       if is_wet_u_ikm1 || is_wet_u_ik # Flux calculated if face is wet and connects to at least one wet u-point
                          w_face = w[i, j, k-1] # w velocity at face k-1 (z_F[k])
                           if w_face > 0 # Downwards flux, upstream is u[i,j,k]
                               u_upstream = is_wet_u_ik ? u[i, j, k] : 0.0 # Should be wet
                           else # Upwards flux, upstream is u[i,j,k-1]
                               u_upstream = is_wet_u_ikm1 ? u[i, j, k-1] : 0.0
                           end
                           flux_n_u_vert = w_face * u_upstream * cc_params.wet_face_area_z[i, j, k-1]
                       end
                end

                # Vertical Divergence = (Flux_into_top - Flux_out_bottom) = (North Flux - South Flux)
                # Divide by the height of the u-cell control volume, which is Δzᵃᵃᶜ[i,j,k].
                Δz_u_cell = grid.Δzᵃᵃᶜ[i,j,k] # Use grid accessor
                if Δz_u_cell > 1e-10
                     # Divergence of vertical flux = (Flux_above - Flux_below) / Δz_u_cell
                     vert_flux_div_u = (flux_n_u_vert - flux_s_u_vert) / Δz_u_cell
                     # Advection adds -(∇ . (uv)). So we subtract the divergence.
                     ∂u∂t.data[i+Hx-1, j+Hy, k+Hz] -= vert_flux_div_u
                end
             end # If hFacW > 0 (using original hFacW)
        end # Loop
        return nothing
    end


    # 5. Vertical Diffusion for u (at Face, Center, Center) - Flux-Based Vertical Viscosity
    # Adds ∂u/∂t = ∂/∂z (Av ∂u/∂z)
    function add_cut_cell_vertical_diffusion_u!(∂u∂t, model)
        grid = model.grid
        u = model.velocities.u
        cc_params = model.parameters.cut_cell_params
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        fill_halo_regions!(u) # Ensure halos are up-to-date for gradient calculation

        # Loop over the interior domain where u-tendency is calculated (F,C,C)
        @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
             # Only compute tendency for wet u-faces (using original hFacW)
             if cc_params.hFacW[i, j, k] > 1e-10

                # Centered difference approximation for ∂/∂z (Av ∂u/∂z) at (i,j,k)
                # Av * [ (u[i,j,k+1] - u[i,j,k]) / Δz_face_k+1 - (u[i,j,k] - u[i,j,k-1]) / Δz_face_k ] / Δz_cell_k
                # Δz_face_k+1 is distance between u[k] and u[k+1] (at face k, index k+1) = grid.Δzᵃᵃᶠ(i,j,k)
                # Δz_face_k is distance between u[k-1] and u[k] (at face k-1, index k) = grid.Δzᵃᵃᶠ(i,j,k-1)
                # Δz_cell_k is height of u[i,j,k] cell = grid.Δzᵃᵃᶜ(i,j,k)

                # Check connectivity using original hFacW
                is_wet_u_ik = cc_params.hFacW[i, j, k] > 1e-10 # u[i,j,k] wet? (Should be true)
                is_wet_u_ikm1 = (k > 1 && cc_params.hFacW[i, j, k-1] > 1e-10) # u[i,j,k-1] wet?
                is_wet_u_ikp1 = (k < grid.Nz && cc_params.hFacW[i, j, k+1] > 1e-10) # u[i,j,k+1] wet?

                # Vertical gradient below u[i,j,k] (between k and k+1)
                grad_below = 0.0
                # Need both u points wet to calculate the gradient between them
                if is_wet_u_ik && is_wet_u_ikp1
                     # Distance between u[i,j,k] and u[i,j,k+1] is Δz at face k+1 (w index k+1), which is grid.Δz_ᵃᵃᶠ(i,j,k)
                     Δz_between = grid.Δz_ᵃᵃᶠ(i,j,k)
                     if Δz_between > 1e-10
                        grad_below = (u[i, j, k+1] - u[i, j, k]) / Δz_between
                     end
                end # If points not connected, gradient is effectively zero

                # Vertical gradient above u[i,j,k] (between k-1 and k)
                grad_above = 0.0
                # Need both u points wet to calculate the gradient between them
                 if is_wet_u_ikm1 && is_wet_u_ik
                     # Distance between u[i,j,k-1] and u[i,j,k] is Δz at face k (w index k), which is grid.Δz_ᵃᵃᶠ(i,j,k-1)
                     Δz_between = grid.Δz_ᵃᵃᶠ(i,j,k-1)
                     if Δz_between > 1e-10
                         grad_above = (u[i, j, k] - u[i, j, k-1]) / Δz_between # Gradient between k-1 and k
                     end
                 end # If points not connected, gradient is effectively zero

                # Vertical Laplacian: Divergence of flux
                # Flux below = -Av * grad_below * Area_face_below (Area_face_z[i,j,k])
                # Flux above = -Av * grad_above * Area_face_above (Area_face_z[i,j,k-1])
                # ∂/∂z(Av ∂u/∂z) approx Av * (grad_above - grad_below) / Δz_u_cell
                # This is simpler than flux divergence here if Av is constant.

                Δz_u_cell = grid.Δzᵃᵃᶜ[i,j,k]
                if Δz_u_cell > 1e-10
                     # Laplacian approx = (grad_above - grad_below) / Δz_u_cell
                    laplacian_z = (grad_above - grad_below) / Δz_u_cell
                    # Add Av * Laplacian to tendency
                    ∂u∂t.data[i+Hx-1, j+Hy, k+Hz] += cc_params.Av * laplacian_z
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
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        Cd = cc_params.bottom_drag_coeff

        fill_halo_regions!(u) # Need halos for u

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
                H_eff = grid.Δzᵃᵃᶜ[i,j,k_bottom] # Use grid accessor

                if H_eff > 1e-10
                    drag_accel = - Cd * abs(u_bottom) * u_bottom / H_eff
                    ∂u∂t.data[i+Hx-1, j+Hy, k_bottom+Hz] += drag_accel
                end
            end
        end
        return nothing
    end


    # Bundle forcings
    u_forcings = (
        pressure_gradient = CustomForcing(add_cut_cell_pressure_gradient_force!, field_dependencies=(:p,)), # PGF depends on pressure
        vertical_advection = CustomForcing(add_cut_cell_vertical_advection_u!, field_dependencies=(:u, :w)), # V-Adv depends on u, w
        vertical_diffusion = CustomForcing(add_cut_cell_vertical_diffusion_u!, field_dependencies=(:u,)),   # V-Diff depends on u
        bottom_drag = CustomForcing(add_cut_cell_bottom_drag!, field_dependencies=(:u,)) # Bottom Drag depends on u
    )

    T_forcings = (
         advection = CustomForcing(add_cut_cell_advection_T!, field_dependencies=(:T, :u, :w)), # Advection depends on T, u, w
         diffusion = CustomForcing(add_cut_cell_diffusion_T!, field_dependencies=(:T,)) # Diffusion depends on T
    )

    S_forcings = (
         advection = CustomForcing(add_cut_cell_advection_T!, field_dependencies=(:S, :u, :w), parameters=nothing), # Reuse T advection function for S
         diffusion = CustomForcing(add_cut_cell_diffusion_S!, field_dependencies=(:S,)) # Custom diffusion for S
    )
    # Note: Reuse T advection for S. The T_forcings advection function takes `(∂T∂t, model)`.
    # We need a generic `add_tracer_advection!(∂c∂t, model)` function, or duplicate the logic.
    # Let's duplicate the logic for now, as the function name should reflect the tracer. Or make it generic.
    # A generic function `add_cut_cell_advection!(∂c∂t, model, tracer_name::Symbol)` is better.
    # Let's refine this.

    # Refined: Make advection generic for any tracer `c`
    function add_cut_cell_advection!(∂c∂t, model)
        grid = model.grid
        c = model.tracers[∂c∂t.name] # Get the tracer field by name

        # Need u and w velocities for advection
        u = model.velocities.u
        w = model.velocities.w # This is the DIAGNOSED w field

        cc_params = model.parameters.cut_cell_params
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        # Ensure required fields have halos filled
        fill_halo_regions!(c)
        fill_halo_regions!(u)
        fill_halo_regions!(w)

        # Loop over the interior grid points where tracer tendency is computed (C,C,C)
        @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
            # Only compute tendency for wet cells (using the merged hFacC mask)
            if cc_params.hFacC[i, j, k] > 1e-10

                # --- Horizontal Fluxes (East - West) of tracer c ---
                # Flux across west face (i) = u[i,j,k] * c_upstream * WetArea_face_x[i,j,k]
                flux_w_c_horiz = 0.0
                if cc_params.hFacW[i, j, k] > 1e-10
                     is_wet_im1k = (i > 1 && cc_params.hFacC[i-1, j, k] > 1e-10)
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10 # Current cell is wet

                     if is_wet_im1k || is_wet_ik
                         u_face = u[i, j, k]
                         if u_face > 0 # Upstream is (i-1, j, k)
                             c_upstream = is_wet_im1k ? c[i-1, j, k] : 0.0
                         else # Upstream is (i, j, k)
                             c_upstream = is_wet_ik ? c[i, j, k] : 0.0
                         end
                         flux_w_c_horiz = u_face * c_upstream * cc_params.wet_face_area_x[i, j, k]
                     end
                end

                # Flux across east face (i+1) = u[i+1,j,k] * c_upstream * WetArea_face_x[i+1,j,k]
                flux_e_c_horiz = 0.0
                if cc_params.hFacW[i+1, j, k] > 1e-10
                     is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                     is_wet_ip1k = (i < grid.Nx && cc_params.hFacC[i+1, j, k] > 1e-10)

                     if is_wet_ik || is_wet_ip1k
                        u_face = u[i+1, j, k]
                        if u_face > 0 # Upstream is (i, j, k)
                            c_upstream = is_wet_ik ? c[i, j, k] : 0.0
                        else # Upstream is (i+1, j, k)
                            c_upstream = is_wet_ip1k ? c[i+1, j, k] : 0.0
                        end
                        flux_e_c_horiz = u_face * c_upstream * cc_params.wet_face_area_x[i+1, j, k]
                     end
                end

                # Horizontal Divergence = (East Flux - West Flux)
                horiz_flux_div_c = flux_e_c_horiz - flux_w_c_horiz

                # --- Vertical Fluxs (North - South) of tracer c ---
                # Flux across face k (at z_F[k+1])
                flux_s_c_vert = 0.0
                 if cc_params.hFacS[i, j, k] > 1e-10
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                      is_wet_ikp1 = (k < grid.Nz && cc_params.hFacC[i, j, k+1] > 1e-10)

                      if is_wet_ik || is_wet_ikp1
                          w_face = w[i, j, k]
                          if w_face > 0 # Downwards flux, upstream is (i,k+1)
                              c_upstream = is_wet_ikp1 ? c[i, j, k+1] : 0.0
                          else # Upwards flux, upstream is (i,k)
                              c_upstream = is_wet_ik ? c[i, j, k] : 0.0
                          end
                          flux_s_c_vert = w_face * c_upstream * cc_params.wet_face_area_z[i, j, k]
                      end
                 end

                # Flux across face k-1 (at z_F[k])
                flux_n_c_vert = 0.0
                if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10
                      is_wet_ikm1 = (k > 1 && cc_params.hFacC[i, j, k-1] > 1e-10)
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10

                       if is_wet_ikm1 || is_wet_ik
                           w_face = w[i, j, k-1]
                           if w_face > 0 # Downwards flux, upstream is (i,k)
                               c_upstream = is_wet_ik ? c[i, j, k] : 0.0
                           else # Upwards flux, upstream is (i,k-1)
                               c_upstream = is_wet_ikm1 ? c[i, j, k-1] : 0.0
                           end
                           flux_n_c_vert = w_face * c_upstream * cc_params.wet_face_area_z[i, j, k-1]
                       end
                end

                # Vertical Divergence = (North Flux - South Flux)
                vert_flux_div_c = flux_n_c_vert - flux_s_c_vert

                # Total Divergence
                total_flux_div_c = horiz_flux_div_c + vert_flux_div_c

                # Tendency = - Total Divergence / Wet Cell Volume (Use merged volume)
                ∂c∂t.data[i+Hx, j+Hy, k+Hz] -= total_flux_div_c * cc_params.recip_wet_cell_volume[i, j, k]

            end # if wet cell (using merged hFacC)
        end # Loop
        return nothing
    end

     # Refined: Make diffusion generic for any tracer `c`
     function add_cut_cell_diffusion!(∂c∂t, model)
         grid = model.grid
         c = model.tracers[∂c∂t.name] # Get the tracer field by name
         cc_params = model.parameters.cut_cell_params
         Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

         fill_halo_regions!(c) # Ensure halos are up-to-date

         @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
             if cc_params.hFacC[i, j, k] > 1e-10

                 # --- Horizontal Fluxes (East - West) of tracer c Gradient ---
                 # Flux across west face (i) = -Kh * grad_x(c) * WetArea_face_x[i,j,k]
                 flux_w_c_horiz = 0.0
                 if cc_params.hFacW[i, j, k] > 1e-10
                      is_wet_im1k = (i > 1 && cc_params.hFacC[i-1, j, k] > 1e-10)
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                      if is_wet_im1k && is_wet_ik
                           Δx_centers = grid.Δx_ᶠᶜᶜ(i,j,k)
                           if Δx_centers > 1e-10
                              grad = (c[i, j, k] - c[i-1, j, k]) / Δx_centers
                              flux_w_c_horiz = -cc_params.Kh * grad * cc_params.wet_face_area_x[i, j, k]
                           end
                      end
                 end

                 # Flux across east face (i+1) = -Kh * grad_x(c) * WetArea_face_x[i+1,j,k]
                 flux_e_c_horiz = 0.0
                 if cc_params.hFacW[i+1, j, k] > 1e-10
                      is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                      is_wet_ip1k = (i < grid.Nx && cc_params.hFacC[i+1, j, k] > 1e-10)
                      if is_wet_ik && is_wet_ip1k
                         Δx_centers = grid.Δx_ᶠᶜᶜ(i+1,j,k)
                         if Δx_centers > 1e-10
                            grad = (c[i+1, j, k] - c[i, j, k]) / Δx_centers
                            flux_e_c_horiz = -cc_params.Kh * grad * cc_params.wet_face_area_x[i+1, j, k]
                         end
                      end
                 end
                 horiz_flux_div_c = flux_e_c_horiz - flux_w_c_horiz

                 # --- Vertical Fluxes (North - South) of tracer c Gradient ---
                 # Flux across face k (at z_F[k+1])
                 flux_s_c_vert = 0.0
                  if cc_params.hFacS[i, j, k] > 1e-10
                       is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                       is_wet_ikp1 = (k < grid.Nz && cc_params.hFacC[i, j, k+1] > 1e-10)

                       if is_wet_ik && is_wet_ikp1
                            Δz_centers = grid.Δz_ᵃᵃᶠ(i,j,k)
                            if Δz_centers > 1e-10
                               grad = (c[i, j, k+1] - c[i, j, k]) / Δz_centers
                               flux_s_c_vert = -cc_params.Kv * grad * cc_params.wet_face_area_z[i, j, k]
                            end
                       end
                  end

                 # Flux across face k-1 (at z_F[k])
                 flux_n_c_vert = 0.0
                 if (k > 1) && cc_params.hFacS[i, j, k-1] > 1e-10
                       is_wet_ikm1 = (k > 1 && cc_params.hFacC[i, j, k-1] > 1e-10)
                       is_wet_ik = cc_params.hFacC[i, j, k] > 1e-10
                       if is_wet_ikm1 && is_wet_ik
                          Δz_centers = grid.Δz_ᵃᵃᶠ(i,j,k-1)
                          if Δz_centers > 1e-10
                              grad = (c[i, j, k] - c[i, j, k-1]) / Δz_centers
                              flux_n_c_vert = -cc_params.Kv * grad * cc_params.wet_face_area_z[i, j, k-1]
                          end
                       end
                 end
                 vert_flux_div_c = flux_n_c_vert - flux_s_c_vert

                 # Total Flux Divergence
                 total_flux_div_c = horiz_flux_div_c + vert_flux_div_c

                 # Tendency = - Total Divergence / Wet Cell Volume (Use merged volume)
                 ∂c∂t.data[i+Hx, j+Hy, k+Hz] -= total_flux_div_c * cc_params.recip_wet_cell_volume[i, j, k]

             end
         end
         return nothing
     end


    # Update the forcings to use the generic functions
    u_forcings = (
        pressure_gradient = CustomForcing(add_cut_cell_pressure_gradient_force!, field_dependencies=(:p,)),
        vertical_advection = CustomForcing(add_cut_cell_vertical_advection_u!, field_dependencies=(:u, :w)),
        vertical_diffusion = CustomForcing(add_cut_cell_vertical_diffusion_u!, field_dependencies=(:u,)),
        bottom_drag = CustomForcing(add_cut_cell_bottom_drag!, field_dependencies=(:u,))
    )

    T_forcings = (
         advection = CustomForcing(add_cut_cell_advection!, field_dependencies=(:T, :u, :w)),
         diffusion = CustomForcing(add_cut_cell_diffusion!, field_dependencies=(:T,))
    )

    S_forcings = (
         advection = CustomForcing(add_cut_cell_advection!, field_dependencies=(:S, :u, :w)),
         diffusion = CustomForcing(add_cut_cell_diffusion!, field_dependencies=(:S,))
    )


    end # module Forcings

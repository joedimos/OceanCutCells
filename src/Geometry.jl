    module Geometry

    using Oceananigans
    using Oceananigans.Units
    using Oceananigans.Grids: xnodes, znodes
    using Oceananigans.Fields: Center, Face, Field

    using ..Parameters: CutCellParameters # Access the struct defined in Parameters.jl

    # --- Bathymetry ---
    function bathymetry_profile(x)
        total_width = 40000.0 # Should get this from parameters ideally
        # For simplicity, hardcoding here or pass parameters to this function if structure is more complex
        base = 500 + 2000 * (x / total_width)
        xc = total_width * 0.6
        w = total_width * 0.08
        ridge = 1500 * exp(-((x - xc)^2) / (2 * w^2))
        depth = base - ridge
        return max(50.0, depth) # Min depth
    end

    function bathymetry_derivative(x)
        total_width = 40000.0 # Hardcoding or pass parameters
        dbase_dx = 2000.0 / total_width
        xc = total_width * 0.6
        w = total_width * 0.08
        # Derivative of exp(-((x - xc)^2) / (2 * w^2)) with respect to x is
        # exp(-((x - xc)^2) / (2 * w^2)) * (-(x - xc) / w^2)
        dridge_dx = 1500.0 * exp(-((x - xc)^2) / (2 * w^2)) * (-(x - xc) / w^2)
        return dbase_dx - dridge_dx
    end

    # --- Initial Cut Cell Geometry Calculation ---
    function compute_initial_geometry(grid, bathymetry_profile)
        Nx, Ny, Nz = grid.Nx, grid.Ny, grid.Nz
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        # Define hFacC, hFacW, and hFacS fields
        hFacC_initial = Field{Center, Center, Center}(grid)
        hFacW_initial = Field{Face, Center, Center}(grid)
        hFacS_initial = Field{Center, Center, Face}(grid)

        # Calculate initial hFac values based on bathymetry
        # x_centers = xnodes(hFacC_initial, Center(), Center(), Center())[:, 1, 1] # xnodes gives full grid
        x_centers_c = xnodes(hFacC_initial, Center(), Center(), Center())[1+Hx:Nx+Hx, 1+Hy:Ny+Hy, 1+Hz:Nz+Hz][:, 1, 1] # Get interior points only
        bathymetry_c = bathymetry_profile.(x_centers_c) # At C locations

        @inbounds for i in 1:Nx, k in 1:Nz
            z_top_face = grid.zᵃᵃᶠ[k] # z_F[k]
            z_bot_face = grid.zᵃᵃᶠ[k+1] # z_F[k+1]
            dz_k = grid.Δzᵃᵃᶜ[i, 1, k] # Δz_C[k] # Using i,1,k indexing for consistency

            bathy_z_at_c = -bathymetry_c[i]

            if z_top_face < bathy_z_at_c # Cell is completely below bathymetry
                hFacC_initial.data[i+Hx, 1+Hy, k+Hz] = 0.0
            elseif z_bot_face >= bathy_z_at_c # Cell is completely above bathymetry
                hFacC_initial.data[i+Hx, 1+Hy, k+Hz] = 1.0
            else # Partial cell
                # hFac = (depth_of_sea_floor - depth_of_cell_bottom_face) / cell_height
                hFacC_initial.data[i+Hx, 1+Hy, k+Hz] = clamp((bathy_z_at_c - z_bot_face) / dz_k, 0.0, 1.0)
            end
        end

        # Simple connectivity masks for faces based on initial hFacC
        # hFacW: 1 if adjacent C cells (i-1,k) and (i,k) are wet, 0 otherwise
        @inbounds for i in 1:Nx+1, k in 1:Nz
            # Use indices i_c for C points
            i_c_left = i - 1
            i_c_right = i

            is_wet_im1k = (i_c_left >= 1 && i_c_left <= Nx) ? hFacC_initial.data[i_c_left+Hx, 1+Hy, k+Hz] > 1e-10 : false
            is_wet_ik = (i_c_right >= 1 && i_c_right <= Nx) ? hFacC_initial.data[i_c_right+Hx, 1+Hy, k+Hz] > 1e-10 : false

            # Solid boundaries at i=1 and i=Nx+1
            if i == 1 || i == Nx + 1
                hFacW_initial.data[i+Hx-1, 1+Hy, k+Hz] = 0.0
            else # Interior face, connects two C cells
                 hFacW_initial.data[i+Hx-1, 1+Hy, k+Hz] = (is_wet_im1k && is_wet_ik) ? 1.0 : 0.0
            end
        end

        # hFacS: 1 if adjacent C cells (i,k) and (i,k+1) are wet, 0 otherwise
        # hFacS[i, j, k] is at Face z index k, which is z_F[k+1]. This face connects cell (i,j,k) and (i,j,k+1).
        @inbounds for i in 1:Nx, k in 1:Nz+1
            # Use indices k_c for C points
            k_c_above = k - 1
            k_c_below = k

            is_wet_ikm1 = (k_c_above >= 1 && k_c_above <= Nz) ? hFacC_initial.data[i+Hx, 1+Hy, k_c_above+Hz] > 1e-10 : false
            is_wet_ik = (k_c_below >= 1 && k_c_below <= Nz) ? hFacC_initial.data[i+Hx, 1+Hy, k_c_below+Hz] > 1e-10 : false

            # Solid bottom boundary at k=Nz+1
            if k == Nz + 1
                hFacS_initial.data[i+Hx, 1+Hy, k+Hz-1] = 0.0
            elseif k == 1 # Open surface (z=0), hFacS[i,j,1] is at z_F[1]
                # This face connects cell (i,j,1) to the boundary above. Assume it's open if the cell below is wet.
                hFacS_initial.data[i+Hx, 1+Hy, k+Hz-1] = (is_wet_ik) ? 1.0 : 0.0 # is_wet_ik here refers to cell (i,j,k) where k=1
            else # Interior face, connects two C cells
                 hFacS_initial.data[i+Hx, 1+Hy, k+Hz-1] = (is_wet_ikm1 && is_wet_ik) ? 1.0 : 0.0
            end
        end


        # Apply hard cut-off
        hFacC_initial[hFacC_initial.data .< 1e-10] .= 0.0
        hFacW_initial[hFacW_initial.data .< 1e-10] .= 0.0
        hFacS_initial[hFacS_initial.data .< 1e-10] .= 0.0

        # Calculate Initial Wet Areas and Volumes (before merging)
        wet_cell_volume_initial = Field{Center, Center, Center}(grid)
        wet_face_area_x_initial = Field{Face, Center, Center}(grid)
        wet_face_area_z_initial = Field{Center, Center, Face}(grid)

        # Wet volume (at C)
        @inbounds for i in 1:Nx, j in 1:Ny, k in 1:Nz
            wet_cell_volume_initial.data[i+Hx, j+Hy, k+Hz] = hFacC_initial[i, j, k] * grid.Δxᶜᵃᵃ[i,j,k] * grid.Δyᵃᶜᵃ[i,j,k] * grid.Δzᵃᵃᶜ[i,j,k] # Use grid Δ volume accessors
        end

        # Wet face area x (at F, west face i)
        @inbounds for i in 1:Nx+1, j in 1:Ny, k in 1:Nz
             # wet_face_area_x is hFacW * dy * dz_at_u (approx Δyᵃᶜᵃ[j] * Δzᵃᵃᶜ[k])
             # But grid areas should be accessed via grid accessor functions like Δyᵃᶜᵃ[i,j,k] and Δzᵃᵃᶜ[i,j,k]
            wet_face_area_x_initial.data[i+Hx-1, j+Hy, k+Hz] = hFacW_initial[i, j, k] * grid.Δyᵃᶜᵃ[i,j,k] * grid.Δzᵃᵃᶜ[i,j,k]
        end

        # Wet face area z (at F, bottom face k)
        @inbounds for i in 1:Nx, j in 1:Ny, k in 1:Nz+1
             # wet_face_area_z is hFacS * dx * dy_at_w (approx Δxᶜᵃᵃ[i] * Δyᵃᶜᵃ[j])
             # But grid areas should be accessed via grid accessor functions like Δxᶜᵃᵃ[i,j,k] and Δyᵃᶜᵃ[i,j,k]
             wet_face_area_z_initial.data[i+Hx, j+Hy, k+Hz-1] = hFacS_initial[i, j, k] * grid.Δxᶜᵃᵃ[i,j,k] * grid.Δyᵃᶜᵃ[i,j,k]
        end


        # Fill halos for initial fields before merging
        fill_halo_regions!(hFacC_initial)
        fill_halo_regions!(hFacW_initial)
        fill_halo_regions!(hFacS_initial)
        fill_halo_regions!(wet_cell_volume_initial)
        fill_halo_regions!(wet_face_area_x_initial)
        fill_halo_regions!(wet_face_area_z_initial)


        return hFacC_initial, hFacW_initial, hFacS_initial,
               wet_cell_volume_initial, wet_face_area_x_initial, wet_face_area_z_initial
    end


    # --- Cell Merging ---
    function merge_cut_cells(grid, initial_hFacC, initial_wet_cell_volume, bathymetry_derivative)
        Nx, Ny, Nz = grid.Nx, grid.Ny, grid.Nz
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        # Copies to modify
        hFacC_merged = deepcopy(initial_hFacC)
        wet_cell_volume_merged = deepcopy(initial_wet_cell_volume)

        # Bathymetry derivative at cell centers in x
        x_centers_c = xnodes(initial_hFacC, Center(), Center(), Center())[1+Hx:Nx+Hx, 1+Hy:Ny+Hy, 1+Hz:Nz+Hz][:, 1, 1] # Get interior points only
        dhdx_c = bathymetry_derivative.(x_centers_c)

        # Mask for cells that will be merged away (source cells)
        source_mask = fill(false, Nx, Ny, Nz)
        # Store the target index for each source cell
        target_idx = fill((-1, -1, -1), Nx, Ny, Nz) # Store (i,j,k) tuple

        # Identify source cells and their targets
        @inbounds for i in 1:Nx, j in 1:Ny, k in 1:Nz
            # Only consider initially wet cells
            if initial_hFacC[i, j, k] > 1e-10
                 regular_cell_volume = grid.Δxᶜᵃᵃ[i,j,k] * grid.Δyᵃᶜᵃ[i,j,k] * grid.Δzᵃᵃᶜ[i,j,k]
                 is_small_cut_cell = (initial_wet_cell_volume[i, j, k] < 0.5 * regular_cell_volume)

                 if is_small_cut_cell
                     # Determine merge direction based on dh/dx (2D logic)
                     grad_x = dhdx_c[i]

                     i_t, j_t, k_t = i, j, k # Default target is self (should not happen if merged)

                     # Prefer vertical merge upwards first? Or based on derivative?
                     # The original code used abs(grad_x) <= 1.0 for vertical. Let's follow that.
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
                     is_target_valid = (i_t >= 1 && i_t <= Nx && j_t >= 1 && j_t <= Ny && k_t >= 1 && k_t <= Nz)
                     is_target_wet = is_target_valid && (initial_hFacC[i_t, j_t, k_t] > 1e-10)

                     if is_target_wet
                         # Mark this cell as a source cell
                         source_mask[i, j, k] = true
                         target_idx[i, j, k] = (i_t, j_t, k_t)
                     else
                          # If target is invalid or dry, this cell cannot be merged into a wet cell.
                          # It remains a small cell in the merged geometry, potentially causing issues.
                          # A warning is appropriate.
                          @warn "Small cell at ($i, $j, $k) could not be merged (invalid or dry target)."
                          source_mask[i, j, k] = false # Do not merge this cell
                     end
                 end
            end # if initial_hFacC > 1e-10
        end

        # Perform the merging of volumes
        # Iterate over source cells and add their volume to the target
        @inbounds for i in 1:Nx, j in 1:Ny, k in 1:Nz
            if source_mask[i, j, k]
                (i_t, j_t, k_t) = target_idx[i, j, k]

                # Add source volume to target volume (using halo-inclusive indices for Field.data)
                wet_cell_volume_merged.data[i_t+Hx, j_t+Hy, k_t+Hz] += initial_wet_cell_volume[i, j, k]

                # Set source volume to zero
                wet_cell_volume_merged.data[i+Hx, j+Hy, k+Hz] = 0.0

                # Mark source hFacC as zero
                hFacC_merged.data[i+Hx, j+Hy, k+Hz] = 0.0

                # Note: Face areas are NOT merged here. This is a simplification carried over.
            end
        end

        num_merged = sum(source_mask)
        if num_merged > 0
             @info "Cell merging complete. Number of merged cells: $(num_merged)"
        else
             @info "Cell merging complete. No cells were merged."
        end


        # Calculate Reciprocal Wet Volume for Merged Grid
        recip_wet_cell_volume_merged = Field{Center, Center, Center}(grid)
        @inbounds for i in 1:Nx, j in 1:Ny, k in 1:Nz
            # Use the modified volume (handle potential floating point issues near zero)
            vol = wet_cell_volume_merged[i,j,k]
            if vol > 1e-15
                recip_wet_cell_volume_merged.data[i+Hx, j+Hy, k+Hz] = 1.0 / vol
            else
                recip_wet_cell_volume_merged.data[i+Hx, j+Hy, k+Hz] = 0.0 # Reciprocal is 0 for dry cells (including merged-away)
            end
        end

        # Ensure halo regions are filled for merged fields
        fill_halo_regions!(hFacC_merged)
        fill_halo_regions!(wet_cell_volume_merged)
        fill_halo_regions!(recip_wet_cell_volume_merged)

        # Return merged hFacC, merged volume, reciprocal volume, and the *initial* face areas
        return hFacC_merged, wet_cell_volume_merged, recip_wet_cell_volume_merged,
               initial_hFacW, initial_hFacS, initial_wet_face_area_x, initial_wet_face_area_z
    end


    end # module Geometry

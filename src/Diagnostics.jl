    module Diagnostics

    using Oceananigans
    using Oceananigans.Fields: Center, Face, Field
    using Oceananigans.BoundaryConditions: fill_halo_regions

    using ..Parameters: CutCellParameters # Access the struct

    # Custom diagnostic to calculate w (at C,C,F) from u (at F,C,C) using continuity
    # ∇ . u = ∂u/∂x + ∂w/∂z = 0
    # Integrate d(wA_z)/dz = -d(uA_x)/dx from bottom up.
    # (wA_z)_at_face_k = (wA_z)_at_face_k+1 + Net_Horiz_Vol_Flux_INTO_cell_k+1
    # Net_Horiz_Vol_Flux_INTO_cell_k+1 = (u*A_x)_W - (u*A_x)_E for cell k+1
    # w_field[i,j,k] corresponds to velocity at face k (z_F[k+1])
    function diagnose_cut_cell_w!(w_field, model)
        grid = model.grid
        u = model.velocities.u
        cc_params = model.parameters.cut_cell_params # Access bundled parameters
        Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz

        fill_halo_regions!(u) # Ensure u halos are up-to-date

        # Temporary field to store vertical volume flux (w*Area_z)
        vertical_volume_flux = Field{Center, Center, Face}(grid)
        # Initialize flux at the bottommost horizontal face (k=Nz+1) to 0.
        # w_field index k=Nz+1 corresponds to face index Nz+1
        w_field.data[i+Hx, j+Hy, grid.Nz+Hz] = 0.0 # w at bottom boundary is 0
        vertical_volume_flux.data[i+Hx, j+Hy, grid.Nz+Hz] = 0.0 # Vertical volume flux at bottom face is 0

        # Integrate upwards from k_cell = grid.Nz down to 1
        @inbounds for i in 1:grid.Nx, j in 1:grid.Ny
            # Initialize vertical flux at the bottom of the water column for this (i,j) pillar
            vertical_volume_flux_at_face_above = 0.0 # This is the flux at face k_cell + 1 initially

            # Loop k_cell from Nz down to 1 (C cell index)
            # The horizontal divergence is for cell k_cell. This divergence drives the flux difference between face k_cell and face k_cell-1.
            # Net horizontal volume flux INTO cell k_cell: (u*A)_W - (u*A)_E
            # (u*A)_W = u[i,j,k_cell] * wet_face_area_x[i,j,k_cell] # at Face i (u index i), at height of cell k_cell
            # (u*A)_E = u[i+1,j,k_cell] * wet_face_area_x[i+1,j,k_cell] # at Face i+1 (u index i+1), at height of cell k_cell

             # Check hFacW for horizontal flux connectivity
             flux_u_w = cc_params.hFacW[i, j, k_cell] > 1e-10 ? u[i, j, k_cell] * cc_params.wet_face_area_x[i, j, k_cell] : 0.0
             flux_u_e = (i <= grid.Nx && cc_params.hFacW[i+1, j, k_cell] > 1e-10) ? u[i+1, j, k_cell] * cc_params.wet_face_area_x[i+1, j, k_cell] : 0.0
             net_horiz_flux_into_cell = flux_u_w - flux_u_e # Net flux INTO cell (i,j,k_cell)

             # Continuity: Sum of fluxes into cell = 0
             # (uA)_W - (uA)_E + (wA)_Top - (wA)_Bottom = 0
             # (wA)_Bottom = (wA)_Top + (uA)_W - (uA)_E
             # (wA) at face k_cell = (wA) at face k_cell - 1 + Net_Horiz_Flux_INTO_cell_k_cell (integrating downwards)

             # Integrate upwards: (wA)_Top = (wA)_Bottom - (uA)_W + (uA)_E
             # (wA) at face k_cell - 1 = (wA) at face k_cell + Net_Horiz_Flux_INTO_cell_k_cell
             # F_{k-1} = F_k + H_k
             # Flux at face k-1 (w_field index k-1) = Flux at face k (w_field index k) + Net horizontal flux into cell k

             # Let vertical_volume_flux[i,j,k] store the flux at FACE k (w_field index k).
             # vertical_volume_flux[i,j,Nz+1] = 0
             # vertical_volume_flux[i,j,k_cell-1] = vertical_volume_flux[i,j,k_cell] + net_horiz_flux_into_cell for cell k_cell

             # This means we sum up the horizontal flux divergence from the bottom.
             # Let's calculate divergence for cell k_cell first.
             horiz_vol_flux_div_cell = flux_u_e - flux_u_w # Divergence = Out - In

             # Vertical volume flux at face k_cell (w_field index k_cell)
             # vvf[i,j,k_cell] = vvf[i,j,k_cell+1] + horiz_vol_flux_div_cell_k_cell
             # This is integrating divergence upwards.

             # vertical_volume_flux.data[i+Hx, j+Hy, k_cell+Hz-1] will store flux at face k_cell (w index k_cell)
             # vertical_volume_flux.data[i+Hx, j+Hy, k_cell+1+Hz-1] stores flux at face k_cell+1 (w index k_cell+1)

             # Flux at face k_cell = Flux at face k_cell + 1 - Horizontal Divergence in cell k_cell + 1? No.
             # Flux at face k = Flux at face k+1 + Net Horizontal Flux Into Cell k+1
             # Flux at face k (w_field index k) = Flux at face k+1 (w_field index k+1) + Net Horiz Flux Into Cell k+1 (C index k+1)

             # Loop k_cell from Nz down to 1 (C index)
             # Flux at face k_cell (w_field index k_cell) = Flux at face k_cell + 1 (w_field index k_cell + 1) + Net Horiz Flux INTO cell k_cell
             k_face_at = k_cell      # Index for w_field and vertical_volume_flux at face k_cell (z_F[k_cell+1])
             k_face_below = k_cell + 1 # Index for w_field and vertical_volume_flux at face k_cell+1 (z_F[k_cell+2])

             # Get flux from face below (already computed in previous iteration, or 0 for k=Nz+1)
             flux_below = vertical_volume_flux.data[i+Hx, j+Hy, k_face_below+Hz-1] # Flux at face k_cell+1

             # Net horizontal volume flux into cell k_cell (C index)
             # This is defined across faces at the same height as cell k_cell.
             flux_u_w = cc_params.hFacW[i, j, k_cell] > 1e-10 ? u[i, j, k_cell] * cc_params.wet_face_area_x[i, j, k_cell] : 0.0
             flux_u_e = cc_params.hFacW[i+1, j, k_cell] > 1e-10 ? u[i+1, j, k_cell] * cc_params.wet_face_area_x[i+1, j, k_cell] : 0.0
             net_horiz_flux_into_cell = flux_u_w - flux_u_e

             # Vertical volume flux at face k_cell (w_field index k_cell)
             # F_k = F_{k+1} + Net_H_flux_into_cell_k
             flux_at_face = flux_below + net_horiz_flux_into_cell
             vertical_volume_flux.data[i+Hx, j+Hy, k_face_at+Hz-1] = flux_at_face

             # Convert flux to velocity W at face k_cell
             if cc_params.wet_face_area_z[i, j, k_face_at] > 1e-15
                  w_field.data[i+Hx, j+Hy, k_face_at+Hz-1] = flux_at_face / cc_params.wet_face_area_z[i, j, k_face_at]
             else
                  w_field.data[i+Hx, j+Hy, k_face_at+Hz-1] = 0.0
             end

             # Ensure W is zero if the face hFacS is zero (explicit mask)
             if cc_params.hFacS[i, j, k_face_at] < 1e-10
                  w_field.data[i+Hx, j+Hy, k_face_at+Hz-1] = 0.0
             end

        end # Loop over k_cell (Nz down to 1)
    end # Loop over i,j

    # Handle top boundary w (at k=1, z_F[1]). This w should be zero for a rigid lid.
    # The integration calculated the flux and velocity at k=1. If z_F[1] is not the surface, this is an internal face.
    # If z_F[1] is the surface (k=1), then the boundary condition w=0 applies here.
    # The current loop calculates w at k=1 based on flux convergence below.
    # If top BC is w=0, need to enforce it. Let's assume w=0 at z=0 (rigid lid).
    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny
        if grid.zᵃᵃᶠ[1] == 0.0 # If the top face is at z=0
             w_field.data[i+Hx, j+Hy, 1+Hz-1] = 0.0 # Set w at k=1 to 0
         end
    end

        # Fill halos for w field after calculation
        fill_halo_regions!(w_field)

        return nothing
    end

    end # module Diagnostics

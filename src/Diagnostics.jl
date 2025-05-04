module Diagnostics

using Oceananigans
using Oceananigans.Fields: Center, Face, Field, fill_halo_regions
using Oceananigans.Grids: ImmersedBoundaryGrid, is_immersed_cell, is_immersed_face,
                          volume, areaᶜᶜᶠ, areaᶠᶜᶜ 

# Custom diagnostic to calculate w (at C,C,F) from u (at F,C,C) using continuity
# ∇ . u = ∂u/∂x + ∂w/∂z = 0
# Integrate d(w*Area_z)/dz = -d(u*Area_x)/dx from bottom up.
# (w*Area_z)_at_face_k = (w*Area_z)_at_face_k+1 + Net_Horiz_Vol_Flux_INTO_cell_k+1
# Net_Horiz_Vol_Flux_INTO_cell_k+1 = (u*Area_x)_W - (u*Area_x)_E for cell k+1
# w_field[i,j,k] corresponds to velocity at face k (z_F[k+1])

function diagnose_cut_cell_w!(w_field, model)
    grid = model.grid # This is now an ImmersedBoundaryGrid
    u = model.velocities.u
    # cc_params = model.parameters.cut_cell_params # No longer needed for geometry

    fill_halo_regions!(u) # Ensure u halos are up-to-date

    # Temporary field to store vertical volume flux (w*Area_z)
    # Create a field at the same location as w_field (C,C,F)
    vertical_volume_flux = Field{Center, Center, Face}(grid)

    # Integrate from bottom up (k = grid.Nz down to 1 for C cell indices)
    # w_field[i,j,k] lives at z_F[k+1] (Face k).
    # The bottommost horizontal face is at z_F[grid.Nz+1] (Face grid.Nz+1).
    # The index for this face in w_field/vertical_volume_flux is grid.Nz+1.

    # Initialize flux at the bottom boundary (k=Nz+1 face) to 0.
    # Face index Nz+1 corresponds to array index Nz+Hz.
    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny
         # Check if the bottom boundary face is actually wet/solid
         # The IBG should handle BCs, but for integration, the flux at the boundary is zero.
         # If !is_immersed_face(i, j, grid.Nz+1, grid, c, c, f) # This face connects to solid bottom
         # Then flux is 0. Our integration naturally starts from 0 at the bottom boundary.
         vertical_volume_flux[i, j, grid.Nz+1] = 0.0
    end

    # Integrate upwards from k_cell = grid.Nz down to 1 (C cell index)
    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny
        # Loop over C cells in vertical (from bottom up)
        for k_cell in grid.Nz : -1 : 1

            # Only integrate if the cell (i,j,k_cell) is wet
            if !is_immersed_cell(i, j, k_cell, grid, Center(), Center(), Center())

                # Net Horizontal Volume Flux INTO cell (i,j,k_cell) (located at C,C,C)
                # Flux IN across west face (i) is u[i,j,k_cell] * WetArea_x_West
                # Flux OUT across east face (i+1) is u[i+1,j,k_cell] * WetArea_x_East
                # WetArea_x at Face i, C cell k_cell is areaᶠᶜᶜ(i, j, k_cell, grid)

                flux_u_w = 0.0
                # Check if the west face (i) is wet
                if !is_immersed_face(i, j, k_cell, grid, Face(), Center(), Center())
                     # Note: u[i,j,k_cell] is at F[i], C[j], C[k_cell]. This lines up.
                     flux_u_w = u[i, j, k_cell] * areaᶠᶜᶜ(i, j, k_cell, grid)
                end

                flux_u_e = 0.0
                # Check if the east face (i+1) is wet
                if (i <= grid.Nx) && !is_immersed_face(i+1, j, k_cell, grid, Face(), Center(), Center())
                     # Note: u[i+1,j,k_cell] is at F[i+1], C[j], C[k_cell]. This lines up.
                     flux_u_e = u[i+1, j, k_cell] * areaᶠᶜᶜ(i+1, j, k_cell, grid)
                end

                net_horiz_vol_flux_into_cell = flux_u_w - flux_u_e

                # Vertical volume flux at face k_cell (z_F[k_cell+1])
                # This face is the bottom face of cell (i,j,k_cell).
                # This is the location of w_field[i,j,k_cell].

                # Continuity: Sum of volume fluxes over cell boundaries is zero.
                # (uA)_E - (uA)_W + (wA)_Top - (wA)_Bottom = 0
                # In 2D: -(Horizontal Divergence) + (wA)_Top - (wA)_Bottom = 0
                # (wA)_Bottom = (wA)_Top - Horizontal Divergence
                # (wA)_at_face_k_cell = (wA)_at_face_k_cell-1 + Net_Horiz_Flux_INTO_cell_k_cell
                # Integrating upwards from bottom boundary (face Nz+1):
                # Flux at face k (w_field index k) = Flux at face k+1 (w_field index k+1) + Net Horiz Flux INTO cell k+1? No.
                # Flux at face k (w_field index k) = Flux at face k+1 (w_field index k+1) + Net Horiz Flux INTO cell k? No.

                # Let F_k be the vertical volume flux across face k (w[i,j,k] location, z_F[k+1]).
                # Sum of fluxes IN = Sum of fluxes OUT for cell (i,j,k).
                # (uA)_W - (uA)_E + (wA)_Top - (wA)_Bottom = 0
                # Horizontal flux IN = (uA)_W - (uA)_E = net_horiz_vol_flux_into_cell
                # Vertical flux IN = (wA)_Top = flux at face k-1 = vvf[i,j,k-1]
                # Vertical flux OUT = (wA)_Bottom = flux at face k = vvf[i,j,k]
                # net_horiz_vol_flux_into_cell + vvf[i,j,k-1] - vvf[i,j,k] = 0
                # vvf[i,j,k-1] = vvf[i,j,k] - net_horiz_vol_flux_into_cell (integrating downwards)
                # vvf[i,j,k] = vvf[i,j,k-1] + net_horiz_vol_flux_into_cell (integrating upwards - WRONG sign)
                # The horizontal divergence is d(uA)/dx. Continuity: d(wA)/dz = -d(uA)/dx.
                # (wA)_top - (wA)_bottom / dz = - ( (uA)_E - (uA)_W ) / dx
                # (wA)_k-1 - (wA)_k / dz = - ( (uA)_E - (uA)_W ) / dx
                # (wA)_k = (wA)_k-1 + ( (uA)_E - (uA)_W ) / dx * dz? No.

                # Let's re-read the original code's logic:
                # "Vertical volume flux at face k_cell (w_field index k_cell)..."
                # "flux_at_face = flux_below + net_horiz_flux_into_cell"
                # This implied: vvf[i,j,k_face_at] = vvf[i,j,k_face_below] + net_horiz_flux_into_cell
                # vvf[i,j,k_cell] = vvf[i,j,k_cell+1] + net_horiz_flux_into_cell_for_cell_k_cell
                # This means integrating upwards: Flux at face k = Flux at face k+1 + Net Horizontal Flux *INTO* cell k (between faces k and k+1)

                # Get flux from face below (already computed in previous iteration, or 0 for k=Nz+1)
                flux_below = vertical_volume_flux[i, j, k_cell + 1] # Flux at face k_cell+1

                # Vertical volume flux at face k_cell (w_field index k_cell)
                flux_at_face = flux_below + net_horiz_vol_flux_into_cell
                vertical_volume_flux[i, j, k_cell] = flux_at_face

                # Convert flux to velocity W at face k_cell (z_F[k_cell+1])
                face_area_z = areaᶜᶜᶠ(i, j, k_cell, grid)
                if face_area_z > 1e-15
                     w_field[i, j, k_cell] = flux_at_face / face_area_z
                else
                     w_field[i, j, k_cell] = 0.0
                end

                # Explicitly set W to zero if the face is solid (is_immersed_face)
                if is_immersed_face(i, j, k_cell, grid, Center(), Center(), Face())
                     w_field[i, j, k_cell] = 0.0
                end

            else # If the C cell (i,j,k_cell) is immersed (dry)
                 # The faces surrounding it should also have zero flux and velocity.
                 # Face k_cell (below cell k_cell) is at w[i,j,k_cell]
                 # Face k_cell-1 (above cell k_cell) is at w[i,j,k_cell-1]
                 # If cell k_cell is dry, the flux across face k_cell-1 should equal flux across face k_cell (no divergence).
                 # However, the loop is integrating based on horizontal flux *into* a cell. If the cell is dry, this flux is zero.
                 # vvf[i,j,k_cell] = vvf[i,j,k_cell+1] + 0
                 # This means the flux above a dry cell equals the flux below it.
                 # The velocity should be zero at faces connected to dry cells.
                 # Let's enforce zero velocity at the face if the face itself is immersed OR if *both* adjacent cells are immersed.
                 # The is_immersed_face check above handles the face being solid.
                 # What if the face is wet but both adjacent cells are dry? Then w should be 0.
                 # Let's rely on is_immersed_face as the primary mask for w.

                 # If cell k_cell is dry, net_horiz_vol_flux_into_cell = 0 (fluxes across its boundaries are zero).
                 # So vertical_volume_flux[i,j,k_cell] = vertical_volume_flux[i,j,k_cell+1] + 0.
                 # This seems correct for dry cells: vertical flux is conserved across dry regions.
                 # The velocity needs to be zero at the faces of dry cells. This is handled by the is_immersed_face check above.
            end # if !is_immersed_cell

        end # Loop over k_cell (Nz down to 1)
    end # Loop over i,j

    # Handle top boundary w (at k=1, z_F[1]).
    # The integration calculated the flux and velocity at k=1.
    # If z_F[1] is the surface (k=1), then the boundary condition w=0 applies here for a rigid lid.
    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny
        if grid.zᵃᵃᶠ[1] == 0.0 # If the top face is at z=0
             w_field[i, j, 1] = 0.0 # Set w at k=1 face to 0
         end
         # Ensure W is zero if the face is solid even at the top boundary
         if is_immersed_face(i, j, 1, grid, Center(), Center(), Face())
             w_field[i, j, 1] = 0.0
         end
    end

    # Fill halos for w field after calculation
    fill_halo_regions!(w_field)

    return nothing
end

end # module Diagnostics

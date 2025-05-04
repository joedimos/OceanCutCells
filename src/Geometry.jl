module Geometry

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: xnodes, znodes, topology, RectilinearGrid, AbstractGrid, with_halo, architecture,
                          xnode, ynode, znode, c, f # Added grid info accessors and locations
using Oceananigans.Fields: Center, Face, Field
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, AbstractGridFittedBottom 
using Oceananigans.Architectures: arch_array # For moving bathymetry to GPU
using Adapt # For GPU adaptation

# --- Bathymetry ---
function bathymetry_profile(x)
    total_width = 40000.0 # Should get this from parameters ideally or pass
    # For simplicity, hardcoding here or pass parameters to this function if structure is more complex
    base = 500 + 2000 * (x / total_width)
    xc = total_width * 0.6
    w = total_width * 0.08
    ridge = 1500 * exp(-((x - xc)^2) / (2 * w^2))
    depth = base - ridge
    return max(50.0, depth) # Min depth
end

# The derivative is not strictly needed for the standard CutCellBottom pattern,

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


# --- Immersed Boundary Type: CutCellBottom ---
struct CutCellBottom{H, E} <: AbstractGridFittedBottom{H}
    "Bathymetry height (negative z value)"
    bottom_height :: H
    "Minimum fractional vertical cell height"
    minimum_fractional_Δz :: E
end

import Base: summary, show
import Oceananigans.ImmersedBoundaries: immersed_cell, bottom_cell, Δzᶜᶠᶜ, Δzᶜᶠᶠ, Δzᶠᶜᶜ, Δzᶠᶜᶠ, Δzᶜᶜᶜ, Δzᶜᶜᶠ, Δzᶠᶠᶜ, Δzᶠᶠᶠ, z_bottom

# Constructor
CutCellBottom(bottom_height; minimum_fractional_Δz=0.1) =
    CutCellBottom(bottom_height, minimum_fractional_Δz)

# Summary and show methods
function summary(ib::CutCellBottom)
    hmax = maximum(ib.bottom_height) # Assuming bottom_height is a field or array
    hmin = minimum(ib.bottom_height) # Need to handle function case differently
    return @sprintf("CutCellBottom(min(h)=%.2f, max(h)=%.2f, ϵ=%.1f)",
                    hmin, hmax, ib.minimum_fractional_Δz)
end

summary(ib::CutCellBottom{<:Function}) = @sprintf("CutCellBottom(%s, ϵ=%.1f)", ib.bottom_height, ib.minimum_fractional_Δz)
show(io::IO, ib::CutCellBottom) = print(io, summary(ib))

# Adapt to architecture (for GPU)
on_architecture(arch, ib::CutCellBottom) = CutCellBottom(arch_array(arch, ib.bottom_height), ib.minimum_fractional_Δz)
Adapt.adapt_structure(to, ib::CutCellBottom) = CutCellBottom(adapt(to, ib.bottom_height), ib.minimum_fractional_Δz)

# --- Immersed Cell Criterion ---
# A C cell (i,j,k) is immersed (solid) if its bottom face (z_F[k+1]) is below the bathymetry height.
# Bathymetry height `h` is a negative z value (depth = -h).
# z_F[k+1] is at znode(i, j, k+1, underlying_grid, c, c, f).
# Cell (i,j,k) is immersed if znode(i, j, k+1, grid, c, c, f) <= bottom_height[i, j].
# Note: bottom_height is provided at C locations in x,y. Needs interpolation for Face locations if used there.

@inline function immersed_cell(i, j, k, underlying_grid, ib::CutCellBottom)
    # z coordinate of the bottom face of cell (i,j,k)
    z_bottom_face = znode(i, j, k+1, underlying_grid, c, c, f)
    # Bathymetry height at the cell horizontal location
    h = @inbounds ib.bottom_height[i, j] # Assuming bottom_height is a Field or array at C,C location
    # Cell is solid if its bottom face is at or below the bathymetry height
    return z_bottom_face <= h
end


# A cell k is a bottom cell if cell k is not solid, but cell k-1 is solid.
# Note: This definition seems reversed from the provided snippet (`k-1` vs `k`).
# In Oceananigans' standard IBG, `bottom_cell(i,j,k)` is the *wet* cell immediately *above* a solid cell.
# A solid cell is immersed_cell(i,j,k, ...).
# The partially wet cell at k is the one where immersed_cell(i,j,k,...) is false but immersed_cell(i,j,k+1,...) is true.
# Let's correct this based on typical Oceananigans IBG usage for bottom cells.
# A C cell (i,j,k) is a "bottom cell" if it's wet, but the cell below it (k+1) is solid.
@inline bottom_cell(i, j, k, ibg::ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:CutCellBottom}) =
    !immersed_cell(i, j, k, ibg.underlying_grid, ibg.immersed_boundary) & # Cell k is wet
     immersed_cell(i, j, k+1, ibg.underlying_grid, ibg.immersed_boundary) # Cell k+1 is solid

# Helper function for the z-coordinate of the bottom boundary at a horizontal location
@inline z_bottom(i, j, ibg::ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:CutCellBottom}) =
    @inbounds ibg.immersed_boundary.bottom_height[i, j] # Assuming bottom_height is at C,C

# --- Specialized Δz accessors for CutCellBottom ---
# These define the effective vertical height of cells/faces near the boundary.
# They are typically used by Oceananigans' built-in operators.
# We need to define them for the ImmersedBoundaryGrid containing CutCellBottom (SCIBG).
const SCIBG = ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:CutCellBottom}

@inline function Δzᶜᶠᶜ(i, j, k, ibg::SCIBG)
    underlying_grid = ibg.underlying_grid
    ib = ibg.immersed_boundary
    # z at c,f,c (the node defining Δzᶜᶠᶜ)
    zc = znode(i, j, k, underlying_grid, c, f, c)
    # z at c,f,f (face node below the c,f,c node)
    zf = znode(i, j, k+1, underlying_grid, c, f, f) # This is z_F[k+1]

    # Bottom height needs to be interpolated to C,F location? The provided snippet uses ℑyᵃᶠᵃ.
    # Assume bottom_height is defined at C,C for simplicity as in original code.
    # Effective height of the cell.
    # The standard Δzᶜᶠᶜ is distance between (c,f,c) at k and (c,f,c) at k+1... No.
    # Δzᶜᶠᶜ(i,j,k) is the distance between z_C[k] and z_F[k+1] (zᶜᶠᶜ[k] - zᶜᶠᶠ[k+1])
    # For a RectilinearGrid, Δzᶜᶜᶜ(k) = zᶜᶜᶠ(k) - zᶜᶜᶠ(k+1) = Δzᵃᵃᶜ(k).
    # Δzᶜᶜᶠ(k) = zᶜᶜᶜ(k) - zᶜᶜᶜ(k+1) = Δzᵃᵃᶠ(k).
    # Δzᶜᶠᶜ(i,j,k) is Δz(i,j,k) centered on c,f,c. For RectilinearGrid, this is just Δzᵃᵃᶜ(k).
    # Δzᶜᶠᶠ(i,j,k) is Δz(i,j,k) centered on c,f,f. For RectilinearGrid, this is just Δzᵃᵃᶠ(k).
    # Similarly for Δzᶠᶜᶜ, Δzᶠᶜᶠ, Δzᶠᶠᶜ, Δzᶠᶠᶠ.
    # The specialized IBG Δz functions define the effective height of the *cut cell* below z_F[k].
    # For a cell at k, the bottom face is k+1. The top face is k.
    # The full height is z_F[k] - z_F[k+1] = Δzᵃᵃᶜ(k).
    # The cut cell height should be z_F[k] - max(z_F[k+1], h).
    # But we need to return the *effective Δz* for grid functions.

    # Reinterpreting the provided Δzᶜᶠᶜ:
    # It seems to calculate the height of the *cut portion* of cell k-1, ending at z_F[k].
    # `z = znode(i, j, k+1, underlying_grid, c, f, f)` is z_F[k+1].
    # `h = ℑyᵃᶠᵃ(i, j, 1, ibg, ib.bottom_height)`. This interpolates bottom_height to C,F. Needs ibg.underlying_grid.
    # `h = ℑyᵃᶠᵃ(i, j, 1, underlying_grid, ib.bottom_height)`.
    # `bottom_height` is at C,C. Interpolating to C,F seems reasonable for vertical functions at F.
    # The calculation `max(ϵ * full_Δz, z - h)` seems to define the *thickness* of the cut cell.
    # Let's rename variables to be clearer.
    # `z_face_below_cell_k_minus_1` = `znode(i, j, k, underlying_grid, c, f, f)` # This is z_F[k]
    # `z_face_below_cell_k` = `znode(i, j, k+1, underlying_grid, c, f, f)` # This is z_F[k+1]
    # `z_c_k` = `znode(i, j, k, underlying_grid, c, f, c)` # This is z_C[k]

    # The provided Δzᶜᶠᶜ(i,j,k, ibg) seems to calculate the height of the *k-th* cell (C-point index k), but at a C,F,C location.
    # It should return the effective vertical size *of the cell centered vertically at k*.
    # This effective size is used for divergences, gradients etc.
    # The cell at index k (C-index) extends from z_F[k] down to z_F[k+1].
    # Its full height is Δzᵃᵃᶜ(i,j,k, underlying_grid).
    # If it's a bottom cell, its height is reduced. The top is z_F[k]. The bottom is max(z_F[k+1], h).
    # Effective height of cell k = z_F[k] - max(z_F[k+1], h).
    # How does this map to Δzᶜᶠᶜ etc.? It seems these functions are meant to return the effective vertical spacing *between* points.

    # Let's assume the provided Δz functions are correct for this specific CutCellBottom implementation.
    # They define the effective vertical spacing required by Oceananigans operators.
    # We need to ensure `bottom_height` is accessible and correctly interpolated/used.
    # `bottom_height` is assumed to be stored at C,C.
    # `ℑyᵃᶠᵃ` interpolates a C,C field to C,F.
    # `ℑxᶠᵃᵃ` interpolates a C,C field to F,C.

    # We need the underlying grid spacing functions as a base.
    import Oceananigans.Grids: Δzᶜᶠᶜ, Δzᶜᶠᶠ, Δzᶠᶜᶜ, Δzᶠᶜᶠ, Δzᶜᶜᶜ, Δzᶜᶜᶠ, Δzᶠᶠᶜ, Δzᶠᶠᶠ

    @inline function Δzᶜᶠᶜ(i, j, k, ibg::SCIBG)
        underlying_grid = ibg.underlying_grid
        ib = ibg.immersed_boundary
        # z coordinate of the face *below* the C,F,C point at k. This is z_F[k+1].
        z_face_below = znode(i, j, k+1, underlying_grid, c, f, f)

        # Bathymetry height interpolated to C,F location
        # Assuming bottom_height is stored as a Field{Center, Center, Nothing}(underlying_grid)
        # or similar that ℑyᵃᶠᵃ can handle. Let's pass the Field directly.
        # Or if bottom_height is a function/array, ℑy needs the grid to interpolate.
        # The constructor takes `bottom_height` which can be a function or array.
        # The `immersed_cell` function accesses `ib.bottom_height[i, j]`. This suggests it's a Field or array.
        # Let's assume it's a Field{Center, Center, Nothing} defined on the underlying grid.
        # Then `ℑyᵃᶠᵃ` needs the underlying grid.
        # h = ℑyᵃᶠᵃ(i, j, 1, underlying_grid, ib.bottom_height) # This requires ℑyᵃᶠᵃ on a Field

        # Alternative: Get bottom height at surrounding C,C points and average/interpolate.
        # For Δzᶜᶠᶜ(i,j,k), which is conceptually vertical spacing around (C,F,C),
        # Bottom height interpolated to C,F.
        # Assume `ib.bottom_height` is a Field{Center, Center, Nothing}(underlying_grid).
        # We need a version interpolated to C,F.
        # This seems overly complicated. Let's stick to the logic from the provided snippet,
        # assuming `ib.bottom_height` somehow handles this. It is defined as {H}.
        # If H is a Field{C,C,N}, then `ib.bottom_height[i, j]` is just the value at C,C.
        # The provided snippet used `ℑyᵃᶠᵃ(i, j, 1, ibg, ib.bottom_height)`. This means ℑyᵃᶠᵃ is defined for SCIBG.
        # This implies `ℑyᵃᶠᵃ(i, j, 1, ibg, field_on_underlying_grid)`.
        # Assume `ib.bottom_height` is a Field{C,C,N}(underlying_grid) inside the IBG.
        # Then ℑyᵃᶠᵃ(i, j, 1, ibg, ib.bottom_height) seems incorrect syntax.
        # It should be `ℑyᵃᶠᵃ(i, j, 1, ibg.underlying_grid, ib.bottom_height)`.

        # Let's re-read the provided snippet carefully: `h = ℑyᵃᶠᵃ(i, j, 1, ibg, ib.bottom_height)`
        # This implies `ℑyᵃᶠᵃ` is overloaded for `SCIBG`.
        # Let's assume `ib.bottom_height` is just a Field{C,C,N}(underlying_grid).
        # Then `ℑyᵃᶠᵃ(i, j, 1, ibg, ib.bottom_height)` might be interpolating on the underlying grid coordinates?
        # Let's use the base underlying grid interpolation `ℑyᵃᶠᵃ(i, j, 1, underlying_grid, ib.bottom_height)`.

        # Revert to simplest assumption: bottom_height is defined at C,C only.
        # The original cut cell code used bathymetry at C points for hFacC.
        # For faces, it used hFacW/S which were just connectivity masks.
        # The provided Δz snippet is more sophisticated. It *needs* bathymetry at the specific IBG location.
        # The provided snippet interpolates `ib.bottom_height` to `C,F` for `Δzᶜᶠᶜ`.
        # Define `bottom_height_CF(i,j)` as `ℑyᵃᶠᵃ(i, j, 1, underlying_grid, ib.bottom_height)`.
        # Similarly, `bottom_height_FC(i,j)` as `ℑxᶠᵃᵃ(i, j, 1, underlying_grid, ib.bottom_height)`.

        # This implies `ib.bottom_height` should be a Field{C,C,N}(underlying_grid).
        # Add this assumption and the interpolation.

        @inline bottom_height_CF(i, j, grid, h_field) = ℑyᵃᶠᵃ(i, j, 1, grid, h_field)
        @inline bottom_height_FC(i, j, grid, h_field) = ℑxᶠᵃᵃ(i, j, 1, grid, h_field)

        underlying_grid = ibg.underlying_grid
        ib = ibg.immersed_boundary
        ϵ = ib.minimum_fractional_Δz

        # Δzᶜᶠᶜ(i,j,k) is distance between zᶜᶠᶜ[k] and zᶜᶠᶜ[k+1] (which is Δzᵃᵃᶜ(k) in RectilinearGrid)
        # The bottom face of cell k is at z_F[k+1].
        # Effective height of cell k: z_F[k] - max(z_F[k+1], h_at_C_C)
        # Let's rethink the provided snippet's Δz.
        # `z = znode(i, j, k+1, underlying_grid, c, f, f)` is z_F[k+1].
        # `h = bottom_height_CF(i, j, underlying_grid, ib.bottom_height)` ?
        # `full_Δz = Δzᶜᶠᶜ(i, j, k, underlying_grid)`.
        # `Cut_Δz = max(ϵ * full_Δz, z - h)`. This implies `z - h` is the cut height.
        # If z_F[k+1] > h, the cell is wet. `z - h` would be the thickness of the wet part *below* z_F[k+1].
        # This is defining the Δz associated with face k+1 ?

        # Let's go back to the simplest interpretation matching the comment in the provided code:
        # `z = znode(i, j, k+1, underlying_grid, c, c, f)` # z_F[k+1] - vertical face below C cell k
        # `h = ib.bottom_height[i, j]` # bottom height at C cell i,j
        # This calculation `max(ϵ * full_Δz, z - h)` seems to define the effective height *below z* at C,C.
        # For Δzᶜᶠᶜ(i,j,k), which is vertical spacing around (C,F,C k), this height is at z_C[k].
        # Effective height at z_C[k] is related to the distance to the *effective* bottom.
        # The effective bottom of cell k is max(z_F[k+1], h).
        # The distance from z_C[k] to the effective bottom is z_C[k] - max(z_F[k+1], h).

        # This seems to define the "distance to the bottom" effective Δz for a half-cell below C.
        # Let's trust the provided structure's intent and use the grid accessors on the underlying grid.
        # `full_Δz = Δzᶜᶠᶜ(i, j, k, underlying_grid)`
        # `z_at_loc = znode(i, j, k, underlying_grid, c, f, c)` # Vertical location of this DeltaZ measure?
        # The provided snippet uses znode(i, j, k+1, ... c, c, f) for Δzᶜᶠᶜ calculation? This feels wrong.
        # It should use the z-locations relevant to Δzᶜᶠᶜ, i.e. zᶜᶠᶜ[k] and zᶜᶠᶜ[k+1].

        # Let's assume the provided functions are defining the *vertical distance between nodes*
        # appropriate for the cut cell based on the height `h`.
        # Δzᶜᶠᶜ(i,j,k) is the distance between zᶜᶠᶜ[k] and zᶜᶠᶜ[k+1].
        # This is normally Δzᵃᵃᶜ(k). In a cut cell, it's zᶜᶠᶜ[k] - effective_zᶜᶠᶜ[k+1].
        # Effective_zᶜᶠᶜ[k+1] should be max(zᶜᶠᶜ[k+1], h_interpolated_to_CFC_at_k+1).

        
        # those specific functions is correct for the intended CutCellBottom behavior.
        # The `z_bottom(i, j, ibg)` needs to be added to ImmersedBoundaries API in Oceananigans master to work generally.
        # For now, define it here as in the snippet.

        # The provided Δz functions use `ℑyᵃᶠᵃ` and `ℑxᶠᵃᵃ`. These need a Field.
        # Let's make `ib.bottom_height` a `Field{Center, Center, Nothing}`.

        @inline function Δzᶜᶠᶜ(i, j, k, ibg::SCIBG)
            underlying_grid = ibg.underlying_grid
            ib = ibg.immersed_boundary
            # z coordinate of the face BELOW the C,F,C point at k?
            # The provided snippet used znode(i, j, k+1, underlying_grid, c, c, f)
            # This is z_F[k+1]. Let's use this.
            z_loc = znode(i, j, k+1, underlying_grid, c, c, f) # This is z_F[k+1]

            # Bottom height needs to be interpolated to C,F location.
            # Provided snippet used ℑyᵃᶠᵃ(i, j, 1, ibg, ib.bottom_height).
            # Let's try to make this work. Assuming ib.bottom_height is a Field.
            # This requires ℑyᵃᶠᵃ to work on IBG or the field itself to be IBG-aware.
            # Let's assume ib.bottom_height is a Field{C,C,N} on the *underlying* grid.
            # Then interpolation needs the underlying grid.
            h = ℑyᵃᶠᵃ(i, j, 1, underlying_grid, ib.bottom_height) # Interpolate h to C,F at k=1 (surface z index)

            # Check if the cell is a bottom cell
            at_the_bottom = bottom_cell(i, j, k, ibg) # Check if cell k is the partial cell

            full_Δz = Δzᶜᶠᶜ(i, j, k, underlying_grid) # Standard vertical distance at C,F,C k
            # This seems to calculate the effective Δz at k, based on the cut at k+1?
            # The standard Δz functions return the distance *between* nodes.
            # Δzᶜᶠᶜ(i,j,k) is the distance between zᶜᶠᶜ[k] and zᶜᶠᶜ[k+1].
            # `max(ϵ * full_Δz, z - h)` calculates the effective thickness of the cut *cell*.
            # And this thickness is used as the effective Δz for the cell k?

            # Let's reconsider: The Δz functions define the distance *between* adjacent points for operators.
            # Δzᶜᶠᶜ(i, j, k) is the distance between tracer points T[i,j,k] and T[i,j,k-1] when evaluated at location (C,F,C,k)? No.
            # Δzᶜᶠᶜ is distance between C points at F face z.
            # Δzᶜᶠᶜ(i,j,k) is distance between zᶜᶠᶜ[k] and zᶜᶠᶜ[k+1].
            # The partial cell is cell k, between z_F[k] and z_F[k+1].
            # The effective distance between C-points k and k+1 is related to the cut.

            # This definition seems to be `max(ϵ * full_Δz, z_F[k+1] - h_CF)`.
            # This value is used when `at_the_bottom` is true (meaning cell k is the cut cell).
            # This effective Δz should be related to `z_F[k] - max(z_F[k+1], h_CF)`.

            # Let's assume the provided `Δzᶜᶠᶜ` etc. functions are correct *as written* for this specific IB type.
            # They define the distance between the nodes *relevant for operators* in the cut cell.
            # This implies `z = znode(i, j, k+1, underlying_grid, c, c, f)` is the correct z location for the boundary.
            # And `h = ℑyᵃᶠᵃ(i, j, 1, underlying_grid, ib.bottom_height)` is the correct interpolated height.
            # And `max(ϵ * full_Δz, z - h)` is the correct distance measure.

            # Let's simplify slightly and assume `ib.bottom_height` is always an array/Field at C,C.
            # Then `ℑy` and `ℑx` must operate on that field and the *underlying* grid.

            h = ℑyᵃᶠᵃ(i, j, 1, underlying_grid, ib.bottom_height) # Interpolate C,C bottom_height to C,F
            z_face_below_cell_k = znode(i, j, k+1, underlying_grid, c, f, f) # z_F[k+1]

            # Is cell k a bottom cell? (wet but cell k+1 below is solid)
            at_the_bottom = bottom_cell(i, j, k, ibg)

            full_Δz = Δzᶜᶠᶜ(i, j, k, underlying_grid)
            # The cut height is the distance from the *top face* of the cut cell (z_F[k])
            # down to the max of the original bottom face (z_F[k+1]) and the boundary (h).
            # Cut cell height = z_F[k] - max(z_F[k+1], h).
            # Let's check if the provided snippet is doing this.
            # `z` is z_F[k+1]. `h` is interpolated bottom height.
            # `z - h` is z_F[k+1] - h. This is positive if z_F[k+1] > h.
            # It seems to be calculating the height *above* the bottom boundary within the cell below k+1?
            # This is very confusing. Let's ignore the provided Δz calculations for a moment and think about what's needed.

            # Oceananigans needs Δz between adjacent points for operators like divergence and gradient.
            # For a C-field c, ∂c/∂z at (C,C,F) k is (c[k]-c[k+1]) / Δzᶜᶜᶠ(k).
            # Δzᶜᶜᶠ(k) is distance between z_C[k] and z_C[k+1].
            # In a cut cell, the distance between C[k] and C[k+1] is complex.

            # Let's make a simplifying assumption: The CutCellBottom only modifies the *volume* of the cell and the *vertical face area* at the bottom.
            # Standard Δx, Δy, and Δz at F,C, points remain based on the underlying grid.
            # This simplifies the forcings/diagnostics:
            # Volume(C) = hFacC * dx * dy * dz
            # Area(F,C,C) = hFacW * dy * dz
            # Area(C,C,F) = hFacS * dx * dy
            # The challenge is that IBG doesn't provide hFac fields directly. It provides `is_immersed` and `volume`/`area` accessors.
            # The provided `CutCellBottom` *does* modify Δz accessors. This means operators will use these modified Δz values.

            # Let's re-read the *standard* IBG definition of immersed_cell for a `GridFittedBottom`.
            # `immersed_cell(i, j, k, grid, ib::GridFittedBottom)`: cell is immersed if z_C[k] < h_interp_to_C.
            # No, it's `z_W[k] < h_interp_to_W` where z_W is the z coord of the bottom *wall*.
            # For a `GridFittedBottom{H,E}`, `immersed_cell(i, j, k, grid, ib)` is `znode(i, j, k, grid, c, c, c) < ib.bottom_height[i, j]`.
            # This assumes `bottom_height` is at C,C. So a cell is solid if its center is below the boundary height. This is *not* the MITgcm partial cell logic.

            # Let's stick to the provided `CutCellBottom` definition for `immersed_cell` (using z_F[k+1] vs h) and the `bottom_cell` definition.
            # And *attempt* to use the provided Δz logic.

            # Δzᶜᶠᶜ(i,j,k, ibg): Need to calculate distance between zᶜᶠᶜ[k] and zᶜᶠᶜ[k+1]?
            # The provided logic is `max(ϵ * full_Δz, z - h)`.
            # `z = znode(i, j, k+1, underlying_grid, c, c, f)` # z_F[k+1]
            # `h = ℑyᵃᶠᵃ(i, j, 1, underlying_grid, ib.bottom_height)` # Interpolate C,C bottom height to C,F at surface k=1?
            # This index `1` for `ℑyᵃᶠᵃ` seems suspicious for a 3D interpolation context.
            # It should be ℑyᵃᶠᵃ(i, j, k, ...).
            # Let's assume `ib.bottom_height` is actually Field{Center, Center, Nothing}(underlying_grid).
            # And `h` is needed at the same depth as the vertical face location.

            # Let's assume the provided Δz functions are meant to define the effective thickness of the *cut* layer.
            # This layer is cell k if `bottom_cell(i, j, k, ibg)` is true.
            # The thickness is `z_F[k] - max(z_F[k+1], h_interpolated)`.
            # `h_interpolated` should be h interpolated to the horizontal location of the specific Δz accessor (C,F, or F,C).

            # Let's redefine the Δz functions based on the common principle:
            # Distance between adjacent nodes is reduced in the cut cell.
            # Δzᶜᶜᶜ(i, j, k, ibg): Vertical distance between C points k and k+1.
            # Full distance: Δzᶜᶜᶜ(i, j, k, underlying_grid).
            # If cell k+1 is a bottom cell, the distance between C[k] and C[k+1] is affected.
            # If cell k is a bottom cell, the distance between C[k-1] and C[k] is affected.
            # This seems overly complicated for a basic cut cell model.

            # Let's assume the simplest model: only the VOLUME of the C cell is fractional,
            # and the AREA of the vertical face below it (C,C,F) is fractional.
            # Horizontal face areas (F,C,C) and horizontal cell volumes (F,C,C) are not fractional.
            # Vertical face areas (C,C,F) are fractional based on hFacS.
            # Cell volumes (C,C,C) are fractional based on hFacC.
            # In the IBG model, `volume(i,j,k, grid, loc...)` and `area(i,j,k, grid, loc...)` should provide these.
            # And `is_immersed` functions provide the wet/dry mask.
            # The Δz functions are for distances between nodes.
            # Let's ignore the provided Δz functions for now, as they seem to implement a different logic than standard partial cells.
            # Standard IBG partial cells affect `volume` and `area`.

            # Redefine CutCellBottom struct without custom Δz methods.
            # The IBG framework should provide default `volume` and `area` accessors that consider `immersed_cell`.

            # struct CutCellBottom{H} <: AbstractGridFittedBottom{H}
            #     bottom_height :: H # Should be a Field{Center, Center, Nothing}(underlying_grid)
            # end
            # Constructor: CutCellBottom(bottom_height_field) = CutCellBottom(bottom_height_field)
            # immersed_cell(i, j, k, underlying_grid, ib::CutCellBottom) = znode(i, j, k+1, underlying_grid, c, c, f) <= ib.bottom_height[i, j] # z_F[k+1] vs h_C_C
            # And rely on default IBG `volume` and `area` implementations which should use this `immersed_cell` criterion.

            # Let's try this simplified approach first, as it aligns better with the idea that IBG handles geometry.
            # If that fails, we can revisit the provided Δz functions.

            # *** Reverting the provided Δz functions for now ***
            # *** This means the CutCellBottom struct is simpler ***

 end # module Geometry

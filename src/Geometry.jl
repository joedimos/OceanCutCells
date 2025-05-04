module Geometry

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: xnodes, znodes, topology, RectilinearGrid, AbstractGrid, with_halo, architecture,
                          xnode, ynode, znode, c, f
using Oceananigans.Fields: Center, Face, Field
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, AbstractGridFittedBottom 
using Oceananigans.Architectures: arch_array
using Adapt
using Oceananigans.Operators: ℑxᶠᵃᵃ, ℑyᵃᶠᵃ

# --- Constants ---
const ϵ = 1e-10  # Small number for floating point comparisons

# --- Bathymetry Functions ---
"""
    bathymetry_profile(x)

Returns the ocean depth (positive downward) at position x.
Implements a sloping bottom with a Gaussian ridge.
"""
function bathymetry_profile(x)
    total_width = 40000.0
    base = 500 + 2000 * (x / total_width)
    xc = total_width * 0.6
    w = total_width * 0.08
    ridge = 1500 * exp(-((x - xc)^2) / (2 * w^2))
    return max(50.0, base - ridge) # Minimum depth of 50m
end

"""
    bathymetry_derivative(x)

Returns the derivative of the bathymetry profile at position x.
Used for calculating bottom slopes.
"""
function bathymetry_derivative(x)
    total_width = 40000.0
    dbase_dx = 2000.0 / total_width
    xc = total_width * 0.6
    w = total_width * 0.08
    dridge_dx = 1500.0 * exp(-((x - xc)^2) / (2 * w^2)) * (-(x - xc) / w^2)
    return dbase_dx - dridge_dx
end

# --- CutCellBottom Type ---
"""
    CutCellBottom{H, E} <: AbstractGridFittedBottom{H}

Immersed boundary type representing a bottom boundary with partial cells.
Contains:
- `bottom_height`: The bathymetry height (negative z value)
- `minimum_fractional_Δz`: Minimum allowed fractional cell height (ϵ)
"""
struct CutCellBottom{H, E} <: AbstractGridFittedBottom{H}
    bottom_height :: H
    minimum_fractional_Δz :: E
end

# Constructor with default minimum fractional height
CutCellBottom(bottom_height; minimum_fractional_Δz=0.1) =
    CutCellBottom(bottom_height, minimum_fractional_Δz)

# Display methods
function summary(ib::CutCellBottom)
    hmax = maximum(ib.bottom_height)
    hmin = minimum(ib.bottom_height)
    return @sprintf("CutCellBottom(min(h)=%.2f, max(h)=%.2f, ϵ=%.1f)",
                    hmin, hmax, ib.minimum_fractional_Δz)
end

summary(ib::CutCellBottom{<:Function}) = @sprintf("CutCellBottom(%s, ϵ=%.1f)", ib.bottom_height, ib.minimum_fractional_Δz)
show(io::IO, ib::CutCellBottom) = print(io, summary(ib))

# GPU/CPU adaptation
on_architecture(arch, ib::CutCellBottom) = CutCellBottom(arch_array(arch, ib.bottom_height), ib.minimum_fractional_Δz)
Adapt.adapt_structure(to, ib::CutCellBottom) = CutCellBottom(adapt(to, ib.bottom_height), ib.minimum_fractional_Δz)

# --- Immersed Boundary Methods ---
"""
    immersed_cell(i, j, k, underlying_grid, ib::CutCellBottom)

Returns true if cell (i,j,k) is immersed (solid), based on whether its bottom face
is below the bathymetry height.
"""
@inline function immersed_cell(i, j, k, underlying_grid, ib::CutCellBottom)
    z_bottom_face = znode(i, j, k+1, underlying_grid, c, c, f)
    h = @inbounds ib.bottom_height[i, j]
    return z_bottom_face <= h
end

"""
    bottom_cell(i, j, k, ibg::ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:CutCellBottom})

Returns true if cell (i,j,k) is the bottom cell (wet cell immediately above solid cells).
"""
@inline function bottom_cell(i, j, k, ibg::ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:CutCellBottom})
    underlying_grid = ibg.underlying_grid
    ib = ibg.immersed_boundary
    return !immersed_cell(i, j, k, underlying_grid, ib) && 
           immersed_cell(i, j, k+1, underlying_grid, ib)
end

"""
    z_bottom(i, j, ibg::ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:CutCellBottom})

Returns the z-coordinate of the bottom boundary at horizontal location (i,j).
"""
@inline z_bottom(i, j, ibg::ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:CutCellBottom}) =
    @inbounds ibg.immersed_boundary.bottom_height[i, j]

# --- Specialized Δz Accessors ---
const SCIBG = ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:Any, <:CutCellBottom}

# Helper functions for interpolation
@inline bottom_height_CF(i, j, grid, h_field) = ℑyᵃᶠᵃ(i, j, grid, h_field)
@inline bottom_height_FC(i, j, grid, h_field) = ℑxᶠᵃᵃ(i, j, grid, h_field)

"""
    Δzᶜᶠᶜ(i, j, k, ibg::SCIBG)

Returns the effective vertical grid spacing at (Center, Face, Center) location,
accounting for partial cells at the bottom boundary.
"""
@inline function Δzᶜᶠᶜ(i, j, k, ibg::SCIBG)
    underlying_grid = ibg.underlying_grid
    ib = ibg.immersed_boundary
    
    z_face = znode(i, j, k+1, underlying_grid, c, c, f)
    h = bottom_height_CF(i, j, underlying_grid, ib.bottom_height)
    full_Δz = Δzᶜᶠᶜ(i, j, k, underlying_grid)
    
    return bottom_cell(i, j, k, ibg) ? max(ib.minimum_fractional_Δz * full_Δz, z_face - h) : full_Δz
end

"""
    Δzᶠᶜᶜ(i, j, k, ibg::SCIBG)

Returns the effective vertical grid spacing at (Face, Center, Center) location,
accounting for partial cells at the bottom boundary.
"""
@inline function Δzᶠᶜᶜ(i, j, k, ibg::SCIBG)
    underlying_grid = ibg.underlying_grid
    ib = ibg.immersed_boundary
    
    z_face = znode(i, j, k+1, underlying_grid, f, c, f)
    h = bottom_height_FC(i, j, underlying_grid, ib.bottom_height)
    full_Δz = Δzᶠᶜᶜ(i, j, k, underlying_grid)
    
    return bottom_cell(i, j, k, ibg) ? max(ib.minimum_fractional_Δz * full_Δz, z_face - h) : full_Δz
end


end # module Geometry
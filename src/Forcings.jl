module Forcings

using Oceananigans
using Oceananigans.Grids: ImmersedBoundaryGrid, is_immersed_cell, is_immersed_face,
                          volume, areaᶜᶜᶠ, areaᶠᶜᶜ,
                          Δx_ᶠᶜᶜ, Δz_ᵃᵃᶠ, Δzᵃᵃᶜ
using Oceananigans.Operators: ∂xᶠᶜᶜ
using Oceananigans.Fields: Center, Face, fill_halo_regions!
using Oceananigans.Forcings: CustomForcing

const ϵ = 1e-10  # Small number for floating point comparisons

"""
    add_cut_cell_pressure_gradient_force!(∂u∂t, model)

Adds pressure gradient force to u-velocity tendency using immersed boundary grid.
Computes ∂u/∂t = -(1/ρ₀) ∂p/∂x at wet faces only.
"""
function add_cut_cell_pressure_gradient_force!(∂u∂t, model)
    grid = model.grid
    p = model.pressure
    inv_rho0 = 1.0 / model.parameters.cut_cell_params.rho0

    # Compute pressure gradient field
    p_gradient = ∂xᶠᶜᶜ(1.0 * p)

    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
        if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
            ∂u∂t[i, j, k] -= inv_rho0 * p_gradient[i, j, k]
        else
            ∂u∂t[i, j, k] = 0.0
        end
    end
    return nothing
end

"""
    add_cut_cell_advection!(∂c∂t, model)

Adds conservative tracer advection using flux-form upwind scheme on immersed boundary grid.
Computes ∂c/∂t = -∇·(u c) at wet cells only.
"""
function add_cut_cell_advection!(∂c∂t, model)
    grid = model.grid
    c = model.tracers[∂c∂t.name]
    u = model.velocities.u
    w = model.velocities.w

    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
        if !is_immersed_cell(i, j, k, grid, Center(), Center(), Center())
            # Horizontal fluxes
            flux_w = flux_e = 0.0
            flux_s = flux_n = 0.0

            # West face
            if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
                u_face = u[i, j, k]
                c_upstream = u_face > 0 ? 
                    (i > 1 && !is_immersed_cell(i-1, j, k, grid) ? c[i-1, j, k] : 0.0) :
                    (!is_immersed_cell(i, j, k, grid) ? c[i, j, k] : 0.0)
                flux_w = u_face * c_upstream * areaᶠᶜᶜ(i, j, k, grid)
            end

            # East face
            if i < grid.Nx && !is_immersed_face(i+1, j, k, grid, Face(), Center(), Center())
                u_face = u[i+1, j, k]
                c_upstream = u_face > 0 ? 
                    (!is_immersed_cell(i, j, k, grid) ? c[i, j, k] : 0.0) :
                    (!is_immersed_cell(i+1, j, k, grid) ? c[i+1, j, k] : 0.0)
                flux_e = u_face * c_upstream * areaᶠᶜᶜ(i+1, j, k, grid)
            end

            # South face
            if !is_immersed_face(i, j, k, grid, Center(), Center(), Face())
                w_face = w[i, j, k]
                c_upstream = w_face > 0 ? 
                    (k < grid.Nz && !is_immersed_cell(i, j, k+1, grid) ? c[i, j, k+1] : 0.0) :
                    (!is_immersed_cell(i, j, k, grid) ? c[i, j, k] : 0.0)
                flux_s = w_face * c_upstream * areaᶜᶜᶠ(i, j, k, grid)
            end

            # North face
            if k > 1 && !is_immersed_face(i, j, k-1, grid, Center(), Center(), Face())
                w_face = w[i, j, k-1]
                c_upstream = w_face > 0 ? 
                    (!is_immersed_cell(i, j, k, grid) ? c[i, j, k] : 0.0) :
                    (!is_immersed_cell(i, j, k-1, grid) ? c[i, j, k-1] : 0.0)
                flux_n = w_face * c_upstream * areaᶜᶜᶠ(i, j, k-1, grid)
            end

            # Compute divergence
            V_cell = volume(i, j, k, grid)
            if V_cell > ϵ
                ∂c∂t[i, j, k] -= ((flux_e - flux_w) + (flux_n - flux_s)) / V_cell
            else
                ∂c∂t[i, j, k] = 0.0
            end
        end
    end
    return nothing
end

"""
    add_cut_cell_diffusion!(∂c∂t, model)

Adds tracer diffusion using flux-form centered differences on immersed boundary grid.
Computes ∂c/∂t = ∇·(K ∇c) at wet cells only.
"""
function add_cut_cell_diffusion!(∂c∂t, model)
    grid = model.grid
    c = model.tracers[∂c∂t.name]
    Kh = model.parameters.cut_cell_params.Kh
    Kv = model.parameters.cut_cell_params.Kv

    @inbounds for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz
        if !is_immersed_cell(i, j, k, grid)
            # Initialize fluxes
            flux_w = flux_e = 0.0
            flux_s = flux_n = 0.0

            # West face
            if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
                im1_ok = i > 1 && !is_immersed_cell(i-1, j, k, grid)
                i_ok = !is_immersed_cell(i, j, k, grid)
                
                if im1_ok && i_ok
                    Δx = Δx_ᶠᶜᶜ(i, j, k, grid)
                    if Δx > ϵ
                        grad = (c[i, j, k] - c[i-1, j, k]) / Δx
                        flux_w = -Kh * grad * areaᶠᶜᶜ(i, j, k, grid)
                    end
                end
            end

            # East face
            if i < grid.Nx && !is_immersed_face(i+1, j, k, grid, Face(), Center(), Center())
                i_ok = !is_immersed_cell(i, j, k, grid)
                ip1_ok = !is_immersed_cell(i+1, j, k, grid)
                
                if i_ok && ip1_ok
                    Δx = Δx_ᶠᶜᶜ(i+1, j, k, grid)
                    if Δx > ϵ
                        grad = (c[i+1, j, k] - c[i, j, k]) / Δx
                        flux_e = -Kh * grad * areaᶠᶜᶜ(i+1, j, k, grid)
                    end
                end
            end

            # South face
            if !is_immersed_face(i, j, k, grid, Center(), Center(), Face())
                k_ok = !is_immersed_cell(i, j, k, grid)
                kp1_ok = k < grid.Nz && !is_immersed_cell(i, j, k+1, grid)
                
                if k_ok && kp1_ok
                    Δz = Δz_ᵃᵃᶠ(i, j, k, grid)
                    if Δz > ϵ
                        grad = (c[i, j, k+1] - c[i, j, k]) / Δz
                        flux_s = -Kv * grad * areaᶜᶜᶠ(i, j, k, grid)
                    end
                end
            end

            # North face
            if k > 1 && !is_immersed_face(i, j, k-1, grid, Center(), Center(), Face())
                km1_ok = !is_immersed_cell(i, j, k-1, grid)
                k_ok = !is_immersed_cell(i, j, k, grid)
                
                if km1_ok && k_ok
                    Δz = Δz_ᵃᵃᶠ(i, j, k-1, grid)
                    if Δz > ϵ
                        grad = (c[i, j, k] - c[i, j, k-1]) / Δz
                        flux_n = -Kv * grad * areaᶜᶜᶠ(i, j, k-1, grid)
                    end
                end
            end

            # Compute divergence
            V_cell = volume(i, j, k, grid)
            if V_cell > ϵ
                ∂c∂t[i, j, k] -= ((flux_e - flux_w) + (flux_n - flux_s)) / V_cell
            else
                ∂c∂t[i, j, k] = 0.0
            end
        end
    end
    return nothing
end

"""
    add_cut_cell_vertical_advection_u!(∂u∂t, model)

Adds vertical advection of u-velocity using flux-form upwind scheme on immersed boundary grid.
Computes ∂u/∂t = -∂(w u)/∂z at wet faces only.
"""
function add_cut_cell_vertical_advection_u!(∂u∂t, model)
    grid = model.grid
    u = model.velocities.u
    w = model.velocities.w

    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
        if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
            # Initialize fluxes
            flux_below = flux_above = 0.0

            # Flux below (at face k)
            if !is_immersed_face(i, j, k, grid, Center(), Center(), Face())
                w_face = w[i, j, k]
                if w_face > 0
                    u_upstream = (k < grid.Nz && !is_immersed_face(i, j, k+1, grid, Face(), Center(), Center())) ? 
                                 u[i, j, k+1] : 0.0
                else
                    u_upstream = !is_immersed_face(i, j, k, grid, Face(), Center(), Center()) ? 
                                u[i, j, k] : 0.0
                end
                flux_below = w_face * u_upstream * areaᶜᶜᶠ(i, j, k, grid)
            end

            # Flux above (at face k-1)
            if k > 1 && !is_immersed_face(i, j, k-1, grid, Center(), Center(), Face())
                w_face = w[i, j, k-1]
                if w_face > 0
                    u_upstream = !is_immersed_face(i, j, k, grid, Face(), Center(), Center()) ? 
                                u[i, j, k] : 0.0
                else
                    u_upstream = !is_immersed_face(i, j, k-1, grid, Face(), Center(), Center()) ? 
                                u[i, j, k-1] : 0.0
                end
                flux_above = w_face * u_upstream * areaᶜᶜᶠ(i, j, k-1, grid)
            end

            V_u_cell = volume(i, j, k, grid, Face(), Center(), Center())
            if V_u_cell > ϵ
                ∂u∂t[i, j, k] -= (flux_below - flux_above) / V_u_cell
            else
                ∂u∂t[i, j, k] = 0.0
            end
        end
    end
    return nothing
end

"""
    add_cut_cell_vertical_diffusion_u!(∂u∂t, model)

Adds vertical diffusion of u-velocity using centered differences on immersed boundary grid.
Computes ∂u/∂t = ∂/∂z(Av ∂u/∂z) at wet faces only.
"""
function add_cut_cell_vertical_diffusion_u!(∂u∂t, model)
    grid = model.grid
    u = model.velocities.u
    Av = model.parameters.cut_cell_params.Av

    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny, k in 1:grid.Nz
        if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
            # Initialize gradients
            grad_below = grad_above = 0.0

            # Gradient below (between k and k+1)
            if !is_immersed_face(i, j, k, grid, Center(), Center(), Face()) &&
               k < grid.Nz && !is_immersed_face(i, j, k+1, grid, Face(), Center(), Center())
                Δz = Δz_ᵃᵃᶠ(i, j, k, grid)
                if Δz > ϵ
                    grad_below = (u[i, j, k+1] - u[i, j, k]) / Δz
                end
            end

            # Gradient above (between k-1 and k)
            if k > 1 && 
               !is_immersed_face(i, j, k-1, grid, Center(), Center(), Face()) &&
               !is_immersed_face(i, j, k-1, grid, Face(), Center(), Center())
                Δz = Δz_ᵃᵃᶠ(i, j, k-1, grid)
                if Δz > ϵ
                    grad_above = (u[i, j, k] - u[i, j, k-1]) / Δz
                end
            end

            Δz_u_cell = Δzᵃᵃᶜ(i, j, k, grid)
            if Δz_u_cell > ϵ
                ∂u∂t[i, j, k] += Av * (grad_above - grad_below) / Δz_u_cell
            end
        else
            ∂u∂t[i, j, k] = 0.0
        end
    end
    return nothing
end

"""
    add_cut_cell_bottom_drag!(∂u∂t, model)

Adds quadratic bottom drag to u-velocity tendency at bottom-most wet faces.
Computes ∂u/∂t = -Cd |u| u / Δz at bottom cells.
"""
function add_cut_cell_bottom_drag!(∂u∂t, model)
    grid = model.grid
    u = model.velocities.u
    Cd = model.parameters.cut_cell_params.bottom_drag_coeff

    @inbounds for i in 1:grid.Nx+1, j in 1:grid.Ny
        # Find bottom-most wet u-face
        k_bottom = -1
        for k in grid.Nz:-1:1
            if !is_immersed_face(i, j, k, grid, Face(), Center(), Center())
                k_bottom = k
                break
            end
        end

        if k_bottom != -1
            u_bottom = u[i, j, k_bottom]
            H_eff = Δzᵃᵃᶜ(i, j, k_bottom, grid)
            
            if H_eff > ϵ
                drag_accel = -Cd * abs(u_bottom) * u_bottom / H_eff
                ∂u∂t[i, j, k_bottom] += drag_accel
            end
        end
    end
    return nothing
end

# Bundle forcings
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
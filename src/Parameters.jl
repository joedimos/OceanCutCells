module Parameters


# Define the struct containing physical parameters
# Removed all geometry fields (hFac, wet_volume, wet_area)
struct CutCellParameters
    # Physical Parameters
    rho0 :: Float64
    g :: Float64
    Kh :: Float64
    Kv :: Float64
    Ah :: Float64
    Av :: Float64
    bottom_drag_coeff :: Float64
    alpha :: Float64
    beta :: Float64
    T0 :: Float64
    S0 :: Float64

    # Parameters needed for model setup (like wind stress)
    wind_stress_tau_x :: Float64 # Added wind stress

    # Grid Parameters (useful to store simulation-specific settings)
    nx :: Int
    nz :: Int
    total_width :: Float64
    total_depth :: Float64
    x_domain :: Tuple{Float64, Float64}
    z_domain :: Tuple{Float64, Float64} # z_domain is informational, grid is built from z_faces
    dt :: Float64 # Time step

end

end # module Parameters

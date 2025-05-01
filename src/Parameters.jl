    module Parameters

    using Oceananigans.Fields: Center, Face, Field

    # Define the struct containing parameters, including cut-cell geometry fields
    struct CutCellParameters
        # Cut-Cell Geometry Fields (These will hold the *merged* or relevant geometry)
        hFacC :: Field{Center, Center, Center}
        hFacW :: Field{Face, Center, Center}
        hFacS :: Field{Center, Center, Face}
        wet_cell_volume   :: Field{Center, Center, Center}
        recip_wet_cell_volume :: Field{Center, Center, Center}
        wet_face_area_x :: Field{Face, Center, Center} # At West face i
        wet_face_area_z :: Field{Center, Center, Face} # At Bottom face k

        # Physical Parameters
        rho0 :: Float64
        g :: Float64
        Kh :: Float64
        Kv :: Float64
        Ah :: Float64
        Av :: Float64
        bottom_drag_coeff :: Float64
        alpha :: Float64 # Added from main script
        beta :: Float64  # Added from main script
        T0 :: Float64   # Added from main script
        S0 :: Float64   # Added from main script

        # Grid Parameters (might be useful to store here or pass separately)
        nx :: Int # Added grid dimensions
        nz :: Int
        total_width :: Float64
        total_depth :: Float64
        x_domain :: Tuple{Float64, Float64}
        z_domain :: Tuple{Float64, Float64}
        dt :: Float64 # Added time step

    end

    end # module Parameters

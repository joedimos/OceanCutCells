    using Test
    using OceanCutCells.Geometry
    using Oceananigans
    using Oceananigans.Fields: Center, Face
    using Oceananigans.Grids: xnodes, znodes

    @testset "Geometry tests" begin
        # Define a simple test grid
        Nx, Nz = 10, 5
        total_width = 100.0
        total_depth = 50.0
        x_domain = (0.0, total_width)
        z_domain = (-total_depth, 0.0)

        grid = VerticallyStretchedRectilinearGrid(
            topology = (Bounded, Flat, Bounded),
            size = (Nx, Nz),
            x = x_domain,
            z_faces = LinRange(-total_depth, 0.0, Nz+1), # Use uniform spacing for simpler testing
            halo = (2, 2)
        )

        # --- Test compute_initial_geometry ---

        # Define a simple test bathymetry: a flat bottom at depth 25m
        function flat_bathymetry(x)
            return 25.0
        end

        hFacC_init, hFacW_init, hFacS_init,
        wet_cell_volume_init, wet_face_area_x_init, wet_face_area_z_init =
            compute_initial_geometry(grid, flat_bathymetry)

        @testset "Initial Geometry (Flat Bottom)" begin
            # Expected hFacC for flat bottom at 25m in 50m depth grid
            # z_faces = [0.0, -10.0, -20.0, -30.0, -40.0, -50.0]
            # z_centers = [-5.0, -15.0, -25.0, -35.0, -45.0]
            # Bathymetry depth: 25m, so -bathymetry_z = -25.0
            # k=1: z_top=0, z_bot=-10. -10 > -25. hFacC = 1.0
            # k=2: z_top=-10, z_bot=-20. -20 > -25. hFacC = 1.0
            # k=3: z_top=-20, z_bot=-30. z_bot < -25 <= z_top. Partial. hFacC = (-25 - (-30)) / (-20 - (-30)) = 5 / 10 = 0.5
            # k=4: z_top=-30, z_bot=-40. -40 < -25. hFacC = 0.0
            # k=5: z_top=-40, z_bot=-50. -50 < -25. hFacC = 0.0

            expected_hFacC = [1.0, 1.0, 0.5, 0.0, 0.0]
            @test all(hFacC_init[i, 1, k] ≈ expected_hFacC[k] for i in 1:Nx, k in 1:Nz)

            # Expected hFacW: 1.0 for all interior faces connecting wet cells
            # Cells (i,1,1), (i,1,2) are wet. Cell (i,1,3) is partial. Cells (i,1,4), (i,1,5) are dry.
            # hFacW[i,1,k] connects (i-1,k) and (i,k)
            # k=1: hFacW = 1.0
            # k=2: hFacW = 1.0
            # k=3: hFacW = 0.0 (connects wet to partial)
            # k=4: hFacW = 0.0
            # k=5: hFacW = 0.0
            expected_hFacW_k = [1.0, 1.0, 0.0, 0.0, 0.0]
            @test all(hFacW_init[i, 1, k] ≈ expected_hFacW_k[k] for i in 2:Nx, k in 1:Nz) # Check interior faces
            @test all(hFacW_init[1, 1, k] ≈ 0.0 for k in 1:Nz) # West boundary
            @test all(hFacW_init[Nx+1, 1, k] ≈ 0.0 for k in 1:Nz) # East boundary

            # Expected hFacS: 1.0 for faces connecting wet cells vertically
            # hFacS[i,j,k] is at Face z index k (z_F[k+1])
            # k=1: Face at z_F[2] (-10). Connects cell 1 and 2. Both wet. hFacS = 1.0
            # k=2: Face at z_F[3] (-20). Connects cell 2 and 3. Cell 2 wet, cell 3 partial. hFacS = 0.0 (Original logic: was 1.0 if *either* was wet. Let's assume connectivity needs both. Rereading the original code: "is_wet_ikm1 && is_wet_ik" -> Needs both wet).
            # k=3: Face at z_F[4] (-30). Connects cell 3 and 4. Cell 3 partial, cell 4 dry. hFacS = 0.0
            # k=4: Face at z_F[5] (-40). Connects cell 4 and 5. Both dry. hFacS = 0.0
            # k=5: Face at z_F[6] (-50). Bottom boundary. hFacS = 0.0
            # k=0 (index 1 in hFacS data): Face at z_F[1] (0). Connects cell 1 to top boundary. Original code: hFacS[i,j,1] = hFacC_initial[i,j,k] > 1e-10 ? 1.0 : 0.0 where k=1 (cell 1). This means hFacS[i,j,1] is 1 if cell 1 is wet. Cell 1 is wet -> hFacS[i,j,1]=1.0.
            expected_hFacS_k = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0] # Indices 1 to Nz+1
            @test all(hFacS_init[i, 1, k] ≈ expected_hFacS_k[k] for i in 1:Nx, k in 1:Nz+1)

            # Check wet volumes/areas are consistent with hFac values and grid spacing
            # With uniform grid: dx=10, dy=1, dz=10 (for C-cells). dx=10, dy=10 (for face-z), dy=1, dz=10 (for face-x)
            # Δxᶜᵃᵃ[i,j,k] = 10, Δyᵃᶜᵃ[i,j,k] = 1, Δzᵃᵃᶜ[i,j,k] = 10 (for interior C)
            # Δyᵃᶜᵃ[i,j,k] * Δzᵃᵃᶜ[i,j,k] = 10 (for Face x area)
            # Δxᶜᵃᵃ[i,j,k] * Δyᵃᶜᵃ[i,j,k] = 10 (for Face z area) -- assuming Face z area calculation is on cell C indices

            # Wet cell volume (at C): hFacC * dx * dy * dz
            expected_wet_volume = expected_hFacC .* (grid.Δxᶜᵃᵃ[1,1,1] * grid.Δyᵃᶜᵃ[1,1,1] * grid.Δzᵃᵃᶜ[1,1,1]) # Volume of regular cell
             @test all(wet_cell_volume_init[i, 1, k] ≈ expected_wet_volume[k] for i in 1:Nx, k in 1:Nz)

            # Wet face area x (at F, west face i): hFacW * dy * dz_at_u
            # dz_at_u seems to be Δzᵃᵃᶜ[i,j,k] in the original code (height of C cell at the same k index).
            # dy_at_u seems to be Δyᵃᶜᵃ[i,j,k].
            # Δyᵃᶜᵃ[i,j,k] * Δzᵃᵃᶜ[i,j,k] = 1 * 10 = 10 for uniform grid
            expected_wet_area_x_k = expected_hFacW_k[1:Nz] .* (grid.Δyᵃᶜᵃ[1,1,1] * grid.Δzᵃᵃᶜ[1,1,1])
             @test all(wet_face_area_x_init[i, 1, k] ≈ expected_wet_area_x_k[k] for i in 2:Nx, k in 1:Nz)

             # Wet face area z (at F, bottom face k): hFacS * dx * dy_at_w
             # dx_at_w seems to be Δxᶜᵃᵃ[i,j,k] in the original code (width of C cell at same i index)
             # dy_at_w seems to be Δyᵃᶜᵃ[i,j,k].
             # Δxᶜᵃᵃ[i,j,k] * Δyᵃᶜᵃ[i,j,k] = 10 * 1 = 10 for uniform grid
            expected_wet_area_z_k = expected_hFacS_k[1:Nz+1] .* (grid.Δxᶜᵃᵃ[1,1,1] * grid.Δyᵃᶜᵃ[1,1,1])
             @test all(wet_face_area_z_init[i, 1, k] ≈ expected_wet_area_z_k[k] for i in 1:Nx, k in 1:Nz+1)
        end

        # --- Test cell merging ---
        # Use the original complex bathymetry to ensure merging happens
        # Need a derivative function that doesn't rely on hardcoded total_width
        function complex_bathymetry_test(x; total_width=100.0)
             base = 500/40000*total_width + (2000/40000*total_width) * (x / total_width)
             xc = total_width * 0.6
             w = total_width * 0.08
             ridge = (1500/40000*total_width) * exp(-((x - xc)^2) / (2 * w^2))
             depth = base - ridge
             return max(50.0/40000*total_width, depth) # Scale min depth
        end

        function complex_bathymetry_derivative_test(x; total_width=100.0)
             dbase_dx = (2000/40000*total_width) / total_width
             xc = total_width * 0.6
             w = total_width * 0.08
             dridge_dx = (1500/40000*total_width) * exp(-((x - xc)^2) / (2 * w^2)) * (-(x - xc) / w^2)
             return dbase_dx - dridge_dx
        end


        # Recompute initial geometry with complex bathymetry on a finer grid
        # Need grid fine enough to generate small cells
        Nx_complex, Nz_complex = 80, 40 # Finer grid
        total_width_complex = 40000.0
        total_depth_complex = 3000.0
        x_domain_complex = (0.0, total_width_complex)
        z_domain_complex = (-total_depth_complex, 0.0)
        z_stretching_complex = 1.3
        z_faces_func_complex(k) = -total_depth_complex * (k/Nz_complex)^z_stretching_complex

         grid_complex = VerticallyStretchedRectilinearGrid(
             topology = (Bounded, Flat, Bounded),
             size = (Nx_complex, Nz_complex),
             x = x_domain_complex,
             z_faces = z_faces_func_complex.(0:Nz_complex),
             halo = (2, 2)
         )

         hFacC_init_comp, hFacW_init_comp, hFacS_init_comp,
         wet_cell_volume_init_comp, wet_face_area_x_init_comp, wet_face_area_z_init_comp =
             compute_initial_geometry(grid_complex, bathymetry_profile)

         hFacC_merged_comp, wet_cell_volume_merged_comp, recip_wet_cell_volume_merged_comp,
         hFacW_geom_comp, hFacS_geom_comp, wet_face_area_x_geom_comp, wet_face_area_z_geom_comp =
             merge_cut_cells(grid_complex, hFacC_init_comp, wet_cell_volume_init_comp, bathymetry_derivative,
                             hFacW_init_comp, hFacS_init_comp, wet_face_area_x_init_comp, wet_face_area_z_init_comp)


        @testset "Cell Merging (Complex Bathymetry)" begin
            # Check that some cells were marked as merged (hFacC_merged == 0 but hFacC_initial > 0)
            num_initial_wet = count(>(1e-10), hFacC_init_comp.data)
            num_merged_zero_hFacC = count(<(1e-10), hFacC_merged_comp.data[grid_complex.Hx+1:end-grid_complex.Hx, 1+grid_complex.Hy:end-grid_complex.Hy, 1+grid_complex.Hz:end-grid_complex.Hz])
            num_initial_zero_hFacC = count(<(1e-10), hFacC_init_comp.data[grid_complex.Hx+1:end-grid_complex.Hx, 1+grid_complex.Hy:end-grid_complex.Hy, 1+grid_complex.Hz:end-grid_complex.Hz])

            @test num_merged_zero_hFacC > num_initial_zero_hFacC # Some cells were newly zeroed out by merging

            # Check that total wet volume is conserved (approx)
            # Sum initial wet volume
            total_volume_initial = sum(wet_cell_volume_init_comp.data)
            # Sum merged wet volume
            total_volume_merged = sum(wet_cell_volume_merged_comp.data)

            # Merging might slightly alter the total volume due to adding small volumes to larger cells.
            # A strict equality check is too sensitive. Check if they are close.
            @test isapprox(total_volume_initial, total_volume_merged, rtol=1e-6) # Check within a relative tolerance

             # Check that dry cells (hFacC_merged < 1e-10) have zero merged volume
            @test all(wet_cell_volume_merged_comp[i,j,k] ≈ 0.0 for i in 1:Nx_complex, j in 1:1, k in 1:Nz_complex if hFacC_merged_comp[i,j,k] < 1e-10)

             # Check that the reciprocal volume is 0 for dry/merged cells and > 0 otherwise
             @test all(recip_wet_cell_volume_merged_comp[i,j,k] ≈ 0.0 for i in 1:Nx_complex, j in 1:1, k in 1:Nz_complex if hFacC_merged_comp[i,j,k] < 1e-10)
             @test all(recip_wet_cell_volume_merged_comp[i,j,k] > 0.0 for i in 1:Nx_complex, j in 1:1, k in 1:Nz_complex if hFacC_merged_comp[i,j,k] >= 1e-10)

             # Check that hFacW, hFacS, face areas are the same as initial (as they weren't merged)
             @test all(hFacW_geom_comp.data .== hFacW_init_comp.data)
             @test all(hFacS_geom_comp.data .== hFacS_init_comp.data)
             @test all(wet_face_area_x_geom_comp.data .== wet_face_area_x_init_comp.data)
             @test all(wet_face_area_z_geom_comp.data .== wet_face_area_z_init_comp.data)
        end

        # @testset "Forcing Tests" begin ... end # Add forcing tests
        # @testset "Diagnostic Tests" begin ... end # Add diagnostic tests
    end

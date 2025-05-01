    using Test
    using OceanCutCells
    using OceanCutCells.Geometry
    using Oceananigans

    @testset "OceanCutCells tests" begin
        include("test_geometry.jl")
        # include("test_forcings.jl") # Add tests for forcings later
        # include("test_model_setup.jl") # Add tests for model setup later
    end

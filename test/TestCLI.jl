function test_parse_command_line()
    push!(ARGS, "--input-file")
    push!(ARGS, "my-input-file.yaml")
    args = Cthonios._parse_command_line()
    @test args["backend"] == "cpu"
    @test args["input-file"] == "my-input-file.yaml"

    push!(ARGS, "--backend")
    push!(ARGS, "gpu")
    push!(ARGS, "--input-file")
    push!(ARGS, "my-input-file.yaml")
    args = Cthonios._parse_command_line()
    @test args["backend"] == "gpu"
    @test args["input-file"] == "my-input-file.yaml"
end

# TODO eventually test for conflicting bcs
function test_parse_dirichlet_bcs()
    settings = Dict(
        Symbol("dirichlet boundary conditions") => [
            Dict(
                Symbol("fields") => ["displ_x", "displ_y"],
                Symbol("function") => "(x, t) -> sin(t)",
                Symbol("sidesets") => ["yminus_sideset", "yplus_sideset"]
            )
        ]
    )
    bcs = Cthonios._parse_dirichclet_bcs(settings)

    @test length(bcs) == 4

    @test bcs[1].sset_name == Symbol("yminus_sideset")
    @test bcs[2].sset_name == Symbol("yplus_sideset")
    @test bcs[3].sset_name == Symbol("yminus_sideset")
    @test bcs[4].sset_name == Symbol("yplus_sideset")

    @test bcs[1].var_name == Symbol("displ_x")
    @test bcs[2].var_name == Symbol("displ_x")
    @test bcs[3].var_name == Symbol("displ_y")
    @test bcs[4].var_name == Symbol("displ_y")

    X = SVector{2, Float64}(0., 0.)
    for bc in bcs
        @test bc.func(X, 0.) ≈ 0.
        @test bc.func(X, π / 2.) ≈ 1.
    end
end

function test_parse_materials()
    settings = Dict(:materials => Dict(
        :mat_1 => Dict(
            :NeoHookean => Dict(
                :density => 1.0,
                Symbol("Young's modulus") => 68.e9,
                Symbol("Poisson's ratio") => 0.3
            ),
            :Gent => Dict(
                :density => 1.,
                Symbol("Young's modulus") => 68.e9,
                Symbol("Poisson's ratio") => 0.3,
                :Jm => 3.
            )
        ),
        :mat_2 => Dict(
            :NeoHookean => Dict(
                :density => 1.,
                Symbol("Young's modulus") => 5.,
                Symbol("Poisson's ratio") => 0.495
            )
        )
    ))
    models, props = Cthonios._parse_materials(settings)

    @test getproperty(models, :mat_1)[:NeoHookean] == NeoHookean
    @test getproperty(models, :mat_1)[:Gent] == Gent
    @test getproperty(models, :mat_2)[:NeoHookean] == NeoHookean

    temp = getproperty(props, :mat_1)[:NeoHookean]
    @test temp[:density] ≈ 1.
    @test temp[Symbol("Young's modulus")] ≈ 68.e9
    @test temp[Symbol("Poisson's ratio")] ≈ 0.3

    temp = getproperty(props, :mat_1)[:Gent]
    @test temp[:density] ≈ 1.
    @test temp[Symbol("Young's modulus")] ≈ 68.e9
    @test temp[Symbol("Poisson's ratio")] ≈ 0.3
    @test temp[:Jm] ≈ 3.

    temp = getproperty(props, :mat_2)[:NeoHookean]
    @test temp[:density] ≈ 1.
    @test temp[Symbol("Young's modulus")] ≈ 5.
    @test temp[Symbol("Poisson's ratio")] ≈ 0.495
end

function test_parse_mesh()
    settings = Dict(
        :mesh => Dict(
            :type => "UnstructuredMesh",
            Symbol("file name") => "my_file.exo"
        )
    )
    type, file_name = Cthonios._parse_mesh(settings)
    @test type == UnstructuredMesh
    @test file_name == "my_file.exo"
end

function test_parse_time()
    settings = Dict(
        :time => Dict(
            Symbol("start time") => 0.,
            Symbol("end time") => 2.,
            Symbol("number of steps") => 15
        )
    )
    time = Cthonios._parse_time(settings)
    @test time.time_current[1] ≈ 0.
    @test time.time_end[1] ≈ 2.
    @test time.Δt[1] ≈ 2. / 15
end

@testset "Parser" begin
    test_parse_command_line()
    test_parse_dirichlet_bcs()
    test_parse_materials()
    test_parse_mesh()
    test_parse_time()
end

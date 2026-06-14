import FiniteElementContainers as FEC
import FiniteElementContainers.AppTools as AT
import FiniteElementContainers.InputFileParser as IP
using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using LinearAlgebra
using SparseArrays

const PHYSICS_LIBRARY = (
    SolidMechanics(PlaneStrain(), Gent()),
    SolidMechanics(PlaneStrain(), NeoHookean()),
    #
    SolidMechanics(ThreeDimensional(), Gent()),
    SolidMechanics(ThreeDimensional(), NeoHookean())
)
const VALID_APP_NAMES = [
    "eigen",
    "explicit",
    "quasistatic"
]

function _setup_physics(sim, fspace)
    materials_settings = IP.get_nested_block(sim.input_settings.parser, "materials")
    physics_settings = IP.get_nested_block(sim.input_settings.parser, "physics")
    material_assignments = physics_settings["material assignment"]::Dict{String, Any}

    # first do a check to make sure all blocks are present in some material
    assigned_blocks = Dict{String, Dict{String, String}}()
    for (mat_name, mat_assignment) in material_assignments
        temp = mat_assignment::Dict{String, Any}
        blocks = IP.get_string_array(temp, "blocks", sim.input_settings.parser.input_style)
        for block in blocks
            if !haskey(assigned_blocks, block)
                model_name = temp["model"]::String
                assigned_blocks[block] = Dict{String, String}("mat_name" => mat_name, "model_name" => model_name)
            else
                @assert false "Repeat material assignment found."
            end
        end
    end

    # check that all blocks are represented
    mesh_blocks = sim.mesh.element_block_names
    physics = AbstractPhysics[]
    properties = Vector{Float64}[]
    initial_states = Vector{Float64}[]
    for block in mesh_blocks
        @assert haskey(assigned_blocks, block)
        mat_name, model_name = assigned_blocks[block]["mat_name"], assigned_blocks[block]["model_name"]
        params = materials_settings[mat_name]::Dict{String, Any}
        params = params[model_name]::Dict{String, Any}

        println(Core.stdout, "Block = ", block)
        if model_name == "Gent"
            temp_physics = PHYSICS_LIBRARY[1]
            temp_props = create_properties(temp_physics, params)
            temp_state = create_initial_state(temp_physics)
        elseif model_name == "NeoHookean"
            temp_physics = PHYSICS_LIBRARY[2]
            temp_props = create_properties(temp_physics, params)
            temp_state = create_initial_state(temp_physics)
        else
            @assert false
        end
        push!(initial_states, temp_state)
        push!(physics, temp_physics)
        push!(properties, temp_props)
        # push!(nstates, length(temp_state))
    end

    # TODO need to actually fill this later
    sizes = FiniteElementContainers.block_quadrature_sizes(fspace)
    nstates = map(length, initial_states)
    state_old = L2Field(undef, Float64, nstates, sizes)
    state_new = L2Field(undef, Float64, nstates, sizes)
    # loop one more time
    for b in axes(mesh_blocks, 1)
        bv_old = FiniteElementContainers.block_view(state_old, b)
        bv_new = FiniteElementContainers.block_view(state_new, b)
        for e in axes(bv_old, 2)
            for q in axes(bv_old, 1)
                temp_state = initial_states[b]
                for s in axes(temp_state, 1)
                    bv_old[s, q, e] = temp_state[s]
                    bv_new[s, q, e] = temp_state[s]
                end
            end
        end
    end
    return physics, properties, state_old, state_new
end

function _cthonios_app_generic_setup(app_name::String, ::Val{D}, ::Val{N}, ARGS::Vector{String}) where {D, N}
    #####################################
    # need to define some types
    #####################################
    SFT = AT.ScalarExpressionFunction{Float64}
    VFT = AT.VectorExpressionFunction{N, Float64}
    SPT = FEC.CSCMatrix()

    ##################################################
    # Setup app
    ##################################################
    app = AT.App{D, N}(app_name)
    AT.add_cli_arg!(app, "--dimension"; is_required = false, short_name = "-d")
    sim = AT.setup(app, ARGS)

    #####################################
    # setup function space
    #####################################
    V = FunctionSpace{true}(sim.mesh, H1Field, Lagrange)
    u = VectorFunction(V, "displ")
    dof = DofManager{false}(u)
    asm = SparseMatrixAssembler{SPT, true, false}(dof)

    # now need to parse materials
    physics, properties, state_old, state_new = _setup_physics(sim, V)
    physics = SolidMechanics{D, 2, 0, PlaneStrain, NeoHookean}[]
    for _ in 1:FiniteElementContainers.num_blocks(V)
        push!(physics, SolidMechanics(PlaneStrain(), NeoHookean()))
    end

    # create objective
    p = FEC.TypeStableParameters{SFT, VFT}(
        sim.mesh, asm,
        physics, properties, state_old, state_new,
        sim.ics, 
        sim.dbcs, sim.nbcs, sim.pbcs, sim.srcs,
        sim.input_settings.time.time
    )
    u = create_unknowns(asm)
    return asm, u, p
end

function _cthonios_app_generic_run(solver, objective, u, p)
    Cthonios.initialize!(objective, u, p)
    # time_end = p.times.time_end
    # n = 2
    # # Cthonios.gradient(objective, u, p)
    # while FiniteElementContainers.current_time(p.times) < time_end - 1e3 * eps(time_end)
    #     Cthonios.step!(solver, objective, u, p)
    #     n = n + 1
    # end
end

function _cthonios_explicit(ARGS::Vector{String}, dimension::Val{D}) where D
    # TODO need to read from input or figure this out
    CFL = 0.1
    asm, u, p = _cthonios_app_generic_setup("explicit", dimension, dimension, ARGS)
    objective = ExplicitDynamicsObjective(asm, CFL)
    solver = ExplicitSolver(objective, u, p)
    _cthonios_app_generic_run(solver, objective, u, p)
end

function _cthonios_quasistatic(ARGS::Vector{String}, dimension::Val{D}) where D
    asm, p = _cthonios_app_generic_setup("quasistatic", dimension, dimension, ARGS)
    objective = QuasiStaticObjective(asm)
    # Cthonios.initialize!(objective)

end

function _get_app_name(ARGS::Vector{String})
    app_name = ARGS[1]
    @assert app_name in VALID_APP_NAMES
    return app_name
end

function _get_dimension_cli_arg(ARGS::Vector{String})
    if "-d" in ARGS
        id = findfirst(x -> x == "-d", ARGS)
    elseif "--dimension" in ARGS
        id = findfirst(x -> x == "--dimension", ARGS)
    else
        @assert false
    end
    dimension = parse(Int, ARGS[id + 1])
    @assert dimension == 2 || dimension == 3
    return dimension, id
end

function _help_message()
    println(
        Core.stdout, 
        "Cthonios cli usage:\n\n",
        "cthonios <app-name>\n",
        "  -d <number of dimensions>\n",
        "  -i <name of input file>\n",
        "  -l <name of file to log to>\n"
    )
end

function @main(ARGS::Vector{String})
    if ARGS[1] == "-h" || ARGS[1] == "--help"
        _help_message()
        Sys.exit()
    end
    app_name = _get_app_name(ARGS)
    dimension, id = _get_dimension_cli_arg(ARGS)
    # purge the app anem and dimension from ARGS
    # deleteat!(ARGS, [1, id, id + 1])
    println(Core.stdout, "App name  = ", app_name)
    println(Core.stdout, "Dimension = ", dimension)
    if dimension == 2
        dim = Val{2}()
    elseif dimension == 3
        dim = Val{3}()
    end

    if app_name == "explicit"
        _cthonios_explicit(ARGS, dim)
    elseif app_name == "quasistatic"
        _cthonios_quasistatic(ARGS, dim)
    end
    return 0
end

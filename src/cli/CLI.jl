function _parse_command_line()::Dict{String, Any}
    settings = ArgParseSettings()
    @add_arg_table! settings begin
      "--backend", "-b"
        arg_type = String
        default = "cpu"
        help = "Backend to use, e.g. cpu, cuda, rocm"
      "--input-file", "-i"
        arg_type = String
        help = "Path to an input file"
        required = true
    end
    return parse_args(settings)
end

function _parse_dirichclet_bcs(settings)
    bc_settings = settings[Symbol("dirichlet boundary conditions")]
    bcs = []
    for bc in bc_settings
        func = @RuntimeGeneratedFunction(Meta.parse(bc[:function]))
        for field in bc[:fields]
            for sset in bc[:sidesets]
                push!(bcs, DirichletBC(field, sset, func))
            end
        end
    end
    return bcs
end

function _parse_input_file(settings)
    input_file = settings["input-file"]
    return YAML.load_file(input_file; dicttype=Dict{Symbol, Any})
end

function _parse_materials(settings)
    material_settings = settings[:materials]
    material_models = Dict{Symbol, Any}()
    material_model_properties = Dict{Symbol, Any}()
    for (name, material) in material_settings
        material_models[name] = Dict{Symbol, Any}()
        material_model_properties[name] = Dict{Symbol, Any}()
        for key in keys(material)
            if key == :name
                continue
            end

            material_models[name][key] = eval(key)
            material_model_properties[name][key] = material[key]
        end
    end
    material_models = NamedTuple{tuple(keys(material_models)...)}(
        tuple(values(material_models)...)
    )
    material_model_properties = NamedTuple{tuple(keys(material_model_properties)...)}(
        tuple(values(material_model_properties)...)
    )
    return material_models, material_model_properties
end

function _parse_mesh(settings)
    mesh_settings = settings[:mesh]
    type = eval(Symbol(mesh_settings[:type]))
    file_name = mesh_settings[Symbol("file name")]
    # return type(file_name)
    return type, file_name
end  

function _parse_neumann_bcs(settings)
    bc_settings = settings[Symbol("neumann boundary conditions")]
    bcs = []
    for bc in bc_settings
        func = @RuntimeGeneratedFunction(Meta.parse(bc[:function]))
        for field in bc[:fields]
            for sset in bc[:sidesets]
                push!(bcs, NeumannBC(field, sset, func))
            end
        end
    end
    return bcs
end

function _parse_nonlinear_solvers(sim_settings)
    nonlinear_solvers = sim_settings[Symbol("nonlinear solvers")]
    return nonlinear_solvers
end

function _parse_physics(settings, material_models, material_model_properties)
    physics_settings = settings[:physics]
    type = eval(Symbol(physics_settings[:type]))
    formulation = eval(Symbol(physics_settings[Symbol("kinematics formulation")]))()
    materials = physics_settings[Symbol("material assignment")]
    physics = Dict{Symbol, Any}()
    props = Dict{Symbol, Any}()

    for (name, material) in pairs(materials)
        model_sym = Symbol(material[:model])
        model_type = material_models[name][model_sym]
        temp_props = Dict{String, Any}()
        for (prop_name, prop) in material_model_properties[name][model_sym]
            temp_props[String(prop_name)] = prop
        end

        for block in material[:blocks]
            block_sym = Symbol(block)
            physics[block_sym] = type(formulation, material_models[name][model_sym]())
            props[block_sym] = temp_props
        end
    end

    physics = NamedTuple(physics)
    props = NamedTuple(props)
    # props = Cthonios.create_properties(physics, props)
    return physics, props
end

function _parse_single_domain_simulation(sim_settings)
    mesh_type, mesh_file = _parse_mesh(sim_settings)
    times = _parse_time(sim_settings)

    material_models, material_model_properties = _parse_materials(sim_settings)
    @info "Materials:"
    for (name, model, props) in zip(keys(material_models), values(material_models), values(material_model_properties))
        @info "  Material: $name"
        for (k, m, p) in zip(keys(model), values(model), values(props))
            @info "    $k:"
            maxlen = maximum(length, string.(keys(p)))
            for (prop_name, prop) in pairs(p)
                padding = " " ^ (maxlen - length(string(prop_name)))
                @info "      $(prop_name)$padding = $prop"
            end
        end
    end
  
    physics, properties = _parse_physics(sim_settings, material_models, material_model_properties)
    @info "Physics:"
    for (name, physics, props) in zip(keys(physics), values(physics), values(properties))
        @info "  Block: $name"
        @info "    Physics: $physics"
        @info "    Properties: $props"
    end
  
    # TODO initial conditions
    dirichlet_bcs = _parse_dirichclet_bcs(sim_settings)
    @info "Dirichlet Boundary Conditions:"
    for (n, bc) in enumerate(dirichlet_bcs)
        # @info bc
        @info "  Dirichlet BC $n"
        @info "    Entity name   = $(bc.sset_name)"
        if isa(bc.func, RuntimeGeneratedFunction)
            @info "    Function      = $(bc.func.body.args[2])"
        else
            @info "    Function      = $(bc.func)"
        end
        @info "    Variable name = $(bc.var_name)"
    end
    dirichlet_bcs = convert(Vector{DirichletBC}, dirichlet_bcs)
  
    neumann_bcs = _parse_neumann_bcs(sim_settings)
    @info "Neumann Boundary Conditions:"
    for bc in neumann_bcs
        @info bc
    end
    neumann_bcs = convert(Vector{NeumannBC}, neumann_bcs)
  
    # TODO add contact pairs

    # TODO solver
    output_file = splitext(mesh_file)[1] * "-output.$(splitext(mesh_file)[2])"
    sim = SingleDomainSimulation(
        mesh_file, output_file, times,
        physics, properties; 
        dirichlet_bcs=dirichlet_bcs
        # neumann_bcs=neumann_bcs
    )
    return sim
end

function _parse_time(sim_settings)
    time_settings = sim_settings[:time]
    start_time = time_settings[Symbol("start time")]
    end_time = time_settings[Symbol("end time")]
    number_of_steps = time_settings[Symbol("number of steps")]
    time = TimeStepper(start_time, end_time, number_of_steps)
    return time
end

# function (@main)(ARGS)
function cthonios_main(ARGS)
    cli_args = _parse_command_line()
    backend = eval(Symbol(cli_args["backend"]))
    input_settings = _parse_input_file(cli_args)
    @info "Backend    = $(cli_args["backend"])"
    @info "Input file = $(cli_args["input-file"])"
    sim_settings = input_settings[:simulation]
    sim_type = sim_settings[:type]
    sim_obj = sim_settings[Symbol("forward problem objective")]
    @info "  Simulation type      = $sim_type"
    @info "  Simulation objective = $sim_obj"

    if sim_type == "SingleDomainSimulation"
        sim = _parse_single_domain_simulation(sim_settings)
        solvers = _parse_nonlinear_solvers(sim_settings)
        @info "Nonlinear solvers:"
        for (k, v) in solvers
            @info "  $k:"
            @info "    $v"
        end

        requested_solver = sim_settings[Symbol("nonlinear solver")]
        solver_settings = solvers[Symbol(requested_solver)]
        solver_type = eval(Symbol(solver_settings[:type]))
        
        kwargs = Dict{Symbol, Any}()
        for (k, v) in solver_settings
            if k == :type
                continue
            end
            kwargs[k] = v
        end
        kwargs = NamedTuple{tuple(keys(kwargs)...)}(tuple(values(kwargs)...))
        solver = x -> solver_type(x; kwargs...)
        Cthonios.run!(sim, QuasiStaticObjectiveCache, solver)
    else
        @assert false "Unsupported simulation type $sim_type"
    end

    return
end
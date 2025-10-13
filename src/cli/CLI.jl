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
    # for material in material_settings
    #     name = Symbol(material[:name])
    #     material_models[name] = eval(Symbol(material[Symbol("material model")]))
    #     material_model_properties[name] = material[Symbol("material model properties")]
    # end
    for (name, material) in material_settings
        # name = Symbol(material[:name])
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
    return material_models, material_model_properties
end

function _parse_mesh(settings)
    mesh_settings = settings[:mesh]
    type = eval(Symbol(mesh_settings[:type]))
    file_name = mesh_settings[Symbol("file name")]
    return type(file_name)
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

function _parse_physics(settings, material_models, material_model_properties)
    physics_settings = settings[:physics]
    type = eval(Symbol(physics_settings[:type]))
    formulation = eval(Symbol(physics_settings[:formulation]))()
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
    props = Cthonios.create_properties(physics, props)
    return physics, props
end

function _parse_single_domain_simulation(sim_settings)
    mesh = _parse_mesh(sim_settings)
    @info "Mesh:"
    # display(mesh)
  
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
  
    physics, props = _parse_physics(sim_settings, material_models, material_model_properties)
    @info "Physics:"
    for (name, physics, props) in zip(keys(physics), values(physics), values(props))
        @info "  Block: $name"
        @info "    Physics: $physics"
        @info "    Properties: $props"
    end
  
    # objective = _parse_objective(sim_settings)
    # @info "Objective: $objective"
  
    # TODO initial conditions

    # TODO time integrator
  
    # TODO initial conditions
    dirichlet_bcs = _parse_dirichclet_bcs(sim_settings)
    @info "Dirichlet Boundary Conditions:"
    for bc in dirichlet_bcs
        @info bc
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
end

function (@main)(ARGS)
    cli_args = _parse_command_line()
    backend = eval(Symbol(cli_args["backend"]))
    input_settings = _parse_input_file(cli_args)
    @info "Backend    = $(cli_args["backend"])"
    @info "Input file = $(cli_args["input-file"])"
    sim_settings = input_settings[:simulation]
    sim_type = sim_settings[:type]
    sim_obj = sim_settings[:objective]
    @info "  Simulation type      = $sim_type"
    @info "  Simulation objective = $sim_obj"

    if sim_type == "SingleDomainSimulation"
        _parse_single_domain_simulation(sim_settings)
    else
        @assert false "Unsupported simulation type $sim_type"
    end

    return
end
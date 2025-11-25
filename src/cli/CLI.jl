include("Parser.jl")

# function (@main)(ARGS)
function cthonios_main()::Cint
    cli_args = _parse_command_line()
    backend = eval(Symbol(cli_args["backend"]))
    input_settings = _parse_input_file(cli_args)
    @info "Backend    = $(cli_args["backend"])"
    @info "Input file = $(cli_args["input-file"])"
    sim_settings = input_settings[:simulation]
    sim_type = sim_settings[:type]
    
    @info "  Simulation type      = $sim_type"
    # @info "  Simulation objective = $sim_obj"

    if sim_type == "SingleDomainSimulation"
        sim = _parse_single_domain_simulation(sim_settings)
        # TODO
        # sim_obj = sim_settings[Symbol("forward problem objective")]
        sim_obj = QuasiStaticObjective()
        obj_cache = setup_cache(sim_obj, sim)

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
        # solver = x -> solver_type(x; kwargs...)
        solver = TrustRegionSolver(obj_cache; kwargs...)
        # Cthonios.run!(sim, QuasiStaticObjectiveCache, solver)
        Cthonios.run!(obj_cache, solver, sim)
    else
        @assert false "Unsupported simulation type $sim_type"
    end

    return 0
end

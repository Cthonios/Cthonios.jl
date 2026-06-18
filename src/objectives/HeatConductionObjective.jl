struct HeatConductionObjective{RT, RV, A} <: AbstractObjective{RT, RV, A}
    assembler::A
    timer::TimerOutput
end

function HeatConductionObjective(assembler, timer = TimerOutput())
    RT = eltype(assembler.constraint_storage)
    RV = typeof(assembler.constraint_storage)
    return HeatConductionObjective{RT, RV, typeof(assembler)}(assembler, timer)
end

function default_output_settings(::HeatConductionObjective)
    return OutputSettings(; displacement = false, temperature = true)
end

function gradient(o::HeatConductionObjective, u, p)
    assemble_vector!(assembler(o), residual!, u, p)
    return residual(assembler(o))
end

function hessian(o::HeatConductionObjective, u, p)
    assemble_stiffness!(assembler(o), stiffness!, u, p)
    return stiffness(assembler(o))
end

function initialize!(::HeatConductionObjective, u, p)
    return nothing
end

function postprocess!(pp::PostProcessor, output_settings, n, objective::HeatConductionObjective, u, p)
    write_times(pp, n, FEC.current_time(p.times))
    if output_settings.temperature
        write_field(pp, n, ("temperature",), p.field)
    end
end

function setup_assembler(::Type{<:HeatConductionObjective}, mesh)
    fspace = FunctionSpace(mesh, H1Field, Lagrange)
    theta = ScalarFunction(fspace, "temperature")
    dof = DofManager(theta)
    assembler = SparseMatrixAssembler(dof; use_inplace_methods = true)
    return assembler
end

function step!(solver, o::HeatConductionObjective, u, p; verbose = false)
    FEC.update_time!(p)
    FEC.update_bc_values!(p, assembler(o))
    _step_begin_banner(o, p; verbose = verbose)
    solve!(solver, o, u, p)

    copyto!(p.field_old, p.field)
    _step_end_banner(o, p; verbose = verbose)
    return nothing
end

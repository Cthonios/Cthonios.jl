struct StructuralOptimization{
    N,
    T <: Number,
    V <: AbstractVector{T},
    M <: AbstractMatrix{T}
} <: AbstractDesignOptimization
    f::V
    dfdϕ::V
    dfdX::H1Field{T, V, N}
    geometry::ESP.Geometry{T, V, M}
    λ::V
    μ::StateVariableField{T, V}
    opt_name::Symbol
    objective_history::V
end

function StructuralOptimization(csm_file::String, opt_name::Symbol)
    geometry = ESP.Geometry(csm_file)
    p = copy(ESP.get_design_parameters(geometry).values)
    ESP.update!(geometry, p)
    p1 = reshape(geometry.dxdps[:, 1], geometry.num_dim, geometry.num_nodes)
    f = zeros(1)
    dfdϕ = zeros(Float64, length(p))
    dfdX = H1Field(zeros(Float64, geometry.num_dim, size(p1, 2)))
    λ = zeros(0)
    μ = StateVariableField(undef, Float64, 0, [(0, 0)])
    # TODO, need to access state vars somehow, initialize as undef maybe?
    # μ = similar(p, )
    # fill!(λ, zero(eltype(λ)))
    objective_history = zeros(0)
    return StructuralOptimization(f, dfdϕ, dfdX, geometry, λ, μ, opt_name, objective_history)
end

function get_constant_parameters(opt::StructuralOptimization)
    return ESP.get_constant_parameters(opt.geometry)
end

function get_design_parameters(opt::StructuralOptimization)
    return ESP.get_design_parameters(opt.geometry)
end

function get_design_sensitivities(opt::StructuralOptimization)
    return opt.geometry.dxdps
end

function optimize!(
    opt::StructuralOptimization, sim_template, solver, qoi;
    maximize::Bool = false,
    kwargs...
)
    p = get_design_parameters(opt)
    x0 = p.values
    # lb = p.lower_bounds
    # ub = p.upper_bounds
    nl_opt = NLopt.Opt(opt.opt_name, length(x0))
    func = (x, g) -> value_and_gradient!(opt, sim_template, solver, qoi, x, g; kwargs...)
    if maximize
        NLopt.max_objective!(nl_opt, func)
    else
        NLopt.min_objective!(nl_opt, func)
    end
    NLopt.lower_bounds!(nl_opt, p.lbs)
    NLopt.upper_bounds!(nl_opt, p.ubs)
    NLopt.xtol_rel!(nl_opt, 1e-8)
    min_f, min_x, ret = NLopt.optimize(nl_opt, x0)
    @show min_f
    @show min_x
    @show ret
end

function _residual!(
    R_u, R, objective, u, p
)
    asm = objective.assembler
    FEC.assemble_vector_enzyme_safe!(R, asm.vector_pattern, asm.dof, residual!, u, p)
    FEC.extract_field_unknowns!(R_u, asm.dof, R)
    return nothing
end

# TODO need to read state variables as well
function value_and_gradient!(opt::StructuralOptimization, sim_template, solver, qoi, x, g; last_step_only::Bool = false)
    ESP.update!(opt.geometry, x)
    sim = run(sim_template, solver)

    objective = Cthonios.get_objective(sim)
    u, p = Cthonios.get_unknowns_and_parameters(sim)

    # reset/resize
    fill!(opt.f, zero(eltype(opt.f)))
    resize!(opt.dfdX, size(p.coords, 2))
    fill!(opt.dfdX, zero(eltype(opt.dfdX)))
    resize!(opt.λ, length(u))
    fill!(opt.λ, zero(eltype(opt.λ)))
    resize!(opt.μ, FEC.block_sizes(p.state_new))
    fill!(opt.μ, zero(eltype(opt.μ)))

    f = zeros(1)
    df = ones(1)
    scratch = Cthonios.create_scratch(qoi, objective)
    dscratch = Cthonios.create_scratch(qoi, objective)
    du = make_zero(u)
    dp = make_zero(p)
    R_u = create_unknowns(objective)
    R = create_field(objective)
    dR_u = make_zero(R_u)
    dR = make_zero(R)

    ExodusDatabase(sim.output_file, "r") do exo
        times = read_times(exo)
        nt = length(times)

        for n in reverse(1:nt)
            # read necessary exodus fields
            displ = H1Field(read_values(exo, NodalVectorVariable, n, "displ"))
            # update time
            p.times.time_current = times[n]
            # TODO need to hook up time step eventually
            # everything assumes constant time step currently

            # update parameters based on current time
            FEC.update_bc_values!(p, objective.assembler)

            # set solution unknown based on read exodus displacement
            FEC.extract_field_unknowns!(u, objective.assembler.dof, displ)

            # assemble stiffness matrix, eventually use operator maybe
            # TODO how can we make this more general for different
            # objective types? Probably hardcoded for quasistatic right now
            assemble_stiffness!(objective.assembler, stiffness!, u, p)
            K = stiffness(objective.assembler)

            # calculate qoi and gradient

            # this bit of enzyme calculate dF/du
            make_zero!(f)
            df[1] = 1.0
            make_zero!(scratch)
            make_zero!(dscratch)
            make_zero!(du)
            make_zero!(dp)
            autodiff(
                Reverse,
                Cthonios.value!,
                Duplicated(f, df),
                Duplicated(scratch, dscratch),
                Const(qoi),
                Const(objective),
                Duplicated(u, du),
                Duplicated(p, dp)
            )

            opt.f .+= f
            opt.dfdX .+= dp.coords

            # TODO replace with Krylov solver maybe?
            opt.λ .= -K \ du

            # residual sensitivities
            make_zero!(R)
            make_zero!(dR)
            make_zero!(R_u)
            make_zero!(dR_u)
            make_zero!(du)
            make_zero!(dp)
            dR_u .= 1.0
            autodiff(
                Reverse,
                _residual!,
                Duplicated(R_u, opt.λ),
                Duplicated(R, dR),
                Const(objective),
                Duplicated(u, du),
                Duplicated(p, dp)
            )
            # display(dp.coords)
            opt.dfdX .+= dp.coords
            # display()

            if last_step_only
                break
            end
        end
    end

    # now calculate dfdphi
    for n in axes(opt.geometry.dxdps, 2)
        opt.dfdϕ[n] = dot(opt.dfdX, opt.geometry.dxdps[:, n])
    end
    # return nothing
    copyto!(g, opt.dfdϕ)
    push!(opt.objective_history, opt.f[1])
    return opt.f[1]
end

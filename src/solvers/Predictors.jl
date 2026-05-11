abstract type AbstractPredictor end

struct TangentPredictor{A, F, V, P} <: AbstractPredictor
    assembler::A
    U::F
    ΔU::F
    ΔUf::V
    R::F
    Rf::V
    dp::P
end

function TangentPredictor(o, p, timer)
    @timeit timer "TangentPredictor - setup" begin
        dof = DofManager(o.assembler.dof.var; use_condensed = true)
        assembler = SparseMatrixAssembler(dof)
        U = create_field(assembler)
        ΔU = create_field(assembler)
        ΔUf = create_unknowns(o)
        R = create_field(assembler)
        Rf = create_unknowns(o)
        dp = deepcopy(p)
        return TangentPredictor(assembler, U, ΔU, ΔUf, R, Rf, dp)
    end
end

function solve!(predictor::TangentPredictor, objective, Uf, p)
    asm = objective.assembler

    # update times in predictor
    # TODO we can probably clean this up a bit
    predictor.dp.times.time_current = p.times.time_current
    predictor.dp.times.Δt = p.times.Δt
    FiniteElementContainers.update_time!(predictor.dp)

    # update bc values
    FiniteElementContainers.update_bc_values!(
        predictor.dp.dirichlet_bcs, p.coords, predictor.dp.times.time_current
    )
    FiniteElementContainers.update_bc_values!(
        predictor.dp.neumann_bcs, predictor.assembler, p.coords, predictor.dp.times.time_current
    )

    predictor.dp.dirichlet_bcs.bc_cache.vals .= p.dirichlet_bcs.bc_cache.vals .- predictor.dp.dirichlet_bcs.bc_cache.vals
    for (bc, dbc) in zip(values(p.neumann_bcs.bc_caches), values(predictor.dp.neumann_bcs.bc_caches))
        dbc.vals .= bc.vals .- dbc.vals
    end

    # form full tangent so we can get Kfc * ΔUc without forming Kfc
    fill!(predictor.U, 0.0)
    update_field_unknowns!(predictor.U, asm.dof, Uf)
    update_field_dirichlet_bcs!(predictor.U, p.dirichlet_bcs)
    assemble_stiffness!(predictor.assembler, objective.objective.hessian_u, predictor.U, p)
    K = stiffness(predictor.assembler)

    # form right hand side for predictor solve
    update_field_dirichlet_bcs!(predictor.ΔU, predictor.dp.dirichlet_bcs)
    mul!(predictor.R.data, K, predictor.ΔU.data)
    # tack on Neumann bc increment
    FiniteElementContainers.assemble_vector_weakly_enforced_bc!(
        predictor.R, predictor.assembler.dof, predictor.ΔU,
        predictor.dp.coords, predictor.dp.neumann_bcs
    )
    FiniteElementContainers.extract_field_unknowns!(predictor.Rf, asm.dof, predictor.R)

    # get stiffness Kff from previous converged step
    assemble_stiffness!(asm, objective.objective.hessian_u, Uf, p)
    Kff = stiffness(asm)

    # solve and update current solution
    predictor.ΔUf .= Kff \ predictor.Rf
    Uf .+= predictor.ΔUf

    return nothing
end

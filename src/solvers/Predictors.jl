abstract type AbstractPredictor end

struct NoPredictor <: AbstractPredictor
end

function NoPredictor(o, u, p, P = I)
    return NoPredictor()
end

function solve!(::NoPredictor, objective, u, p, P)
    return nothing
end

struct TangentPredictor{A, F, V, P, S} <: AbstractPredictor
    assembler::A
    U::F
    ΔU::F
    ΔUf::V
    R::F
    Rf::V
    dp::P
    linear_solver::S
    timer::TimerOutput
end

function TangentPredictor(o, u, p, solver_type = Val{:cg}(), timer = TimerOutput())
    @timeit timer "TangentPredictor - setup" begin
        dof = DofManager{true}(o.assembler.dof.var)
        asm = SparseMatrixAssembler(
            dof;
            matrix_free = true,
            use_inplace_methods = FiniteElementContainers._use_inplace_methods(o.assembler)
        )
        U = create_field(asm)
        ΔU = create_field(asm)
        ΔUf = create_unknowns(o)
        R = create_field(asm)
        Rf = create_unknowns(o)
        dp = deepcopy(p)

        # setup krylov solver workspace
        Kff = stiffness(o.assembler)
        linear_solver = KrylovSolver(solver_type, Kff, Rf; timer = timer)

        return TangentPredictor(asm, U, ΔU, ΔUf, R, Rf, dp, linear_solver, timer)
    end
end

function solve!(predictor::TangentPredictor, objective, Uf, p, P)
    @timeit predictor.timer "TangentPredictor - solve!" begin
        asm = objective.assembler

        # update times in predictor
        # TODO we can probably clean this up a bit
        predictor.dp.times.time_current = p.times.time_current
        predictor.dp.times.Δt = p.times.Δt
        FiniteElementContainers.update_time!(predictor.dp)

        # update bc values
        @timeit predictor.timer "update bc values" begin
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
        end

        # form full tangent so we can get Kfc * ΔUc without forming Kfc
        @timeit predictor.timer "update field" begin
            fill!(predictor.U, 0.0)
            update_field_unknowns!(predictor.U, asm.dof, Uf)
            update_field_dirichlet_bcs!(predictor.U, p.dirichlet_bcs)
            # form right hand side for predictor solve
            update_field_dirichlet_bcs!(predictor.ΔU, predictor.dp.dirichlet_bcs)
        end

        @timeit predictor.timer "assembly" begin
            assemble_matrix_action!(
                predictor.R, predictor.assembler.vector_pattern, predictor.assembler.dof, 
                stiffness_action!, predictor.U, predictor.ΔU, predictor.dp;
                use_inplace_methods = FiniteElementContainers._use_inplace_methods(predictor.assembler)
            )

            # tack on Neumann bc increment
            FiniteElementContainers.assemble_vector_weakly_enforced_bc!(
                predictor.R, predictor.assembler.dof, predictor.ΔU,
                predictor.dp.coords, predictor.dp.neumann_bcs
            )
            # TODO need to add body forces / point loads
            FiniteElementContainers.extract_field_unknowns!(predictor.Rf, asm.dof, predictor.R)

            # get stiffness Kff from previous converged step
            assemble_stiffness!(asm, stiffness!, Uf, p)
            Kff = stiffness(asm)
        end

        # solve and update current solution
        ΔUf = solve!(predictor.linear_solver, Kff, predictor.Rf, P)
        Uf .+= ΔUf
    end
    return nothing
end

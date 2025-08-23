struct CauchyPoint{C, T}
    cauchy_point::C
    timer::T
end

function calculate!(
    cauchy_point::CauchyPoint, 
    objective_cache,
    Uu,
    p,
    g, 
    # hess_vec_func,
    mult_by_approx_hessian,
    tr_size
)
    @timeit cauchy_point.timer "Cauchy Point" begin
        cp = cauchy_point.cauchy_point
        Kg = hvp(objective_cache, Uu, p, g)
        # gKg = dot(g, hess_vec_func(g))
        gKg = dot(g, Kg)
        if gKg > 0
            alpha = -dot(g, g) / gKg
            # cauchyPoint = alpha * g
            copyto!(cp, alpha * g)
            cp_norm_sq = dot(cp, mult_by_approx_hessian(cp))
        else
            # cauchyPoint =  -g * (tr_size / sqrt.(dot(g, mult_by_approx_hessian(g))))
            copyto!(cp, -g * (tr_size / sqrt.(dot(g, mult_by_approx_hessian(g)))))
            cp_norm_sq = tr_size * tr_size
            @warn "negative curvature unpreconditioned cauchy point direction found."
        end

        return cp_norm_sq
    end
end

struct DogLegStep{T}
    timer::T
end

function solve!(::DogLegStep, cp, newton_p, tr_size, mat_mul)
    cc = dot(cp, mat_mul(cp))
    nn = dot(newtonP, mat_mul(newton_p))
    tt = tr_size * tr_size
      
    if cc >= tt #return cauchy point if it extends outside the tr
        #print('cp on boundary')
        return cp * sqrt(tt/cc)
    end
  
    if cc > nn # return cauchy point?  seems the preconditioner was not accurate?
        @warn "cp outside newton, preconditioner likely inaccurate"
        return cp
    end
  
    if nn > tt # on the dogleg (we have nn >= cc, and tt >= cc)
        return preconditioned_project_to_boundary(
            cp,
            newtonP-cp,
            trSize,
            cc,
            mat_mul
        )
    end
    return newtonP  
end

struct TrustRegionModelProblem{V, T}
    cp::V
    d::V
    Pr::V
    z::V
    timer::T
end

function minimize!(
    model_problem::TrustRegionModelProblem,
    x, r, P, hess_vec_func, tr_size, settings
)
    @timeit model_problem.timer "Model Problem" begin
        cp = model_problem.cp
        d = model_problem.d
        Pr = model_problem.Pr
        z = model_problem.z

        fill!(z, zero(eltype(x))) # z = 0. * x
        zz = 0.

        cg_inexact_rel_tol = settings.cg_inexact_solve_ratio
        cg_tol_sq = max(settings.cg_tol^2, cg_inexact_rel_tol * cg_inexact_rel_tol * dot(r, r))
        if dot(r, r) < cg_tol_sq
            return z, z, interiorString, 0
            return nothing
        end

        ldiv!(Pr, P, r)
        copyto!(d, -Pr)
        copyto!(cp, Pr)
        rPr = dot(r, Pr)

        zz = 0.0
        zd = 0.0

        if settings.use_preconditioned_inner_product_for_cg
            dd = rPr
            cg_inner_products = cg_inner_products_preconditioned
        else
            dd = dot(d, d)
            cg_inner_products = cg_inner_products_unpreconditioned
        end

        # begin cg iters
        for i in 1:settings.max_cg_iters
            curvature = dot(d, hess_vec_func(d))
            alpha = rPr / curvature
                
            zNp1 = z + alpha*d
            zzNp1 = update_step_length_squared(alpha, zz, zd, dd)
        
            if curvature <= 0
                zOut = project_to_boundary_with_coefs(z, d, tr_size,
                                                      zz, zd, dd)
                return zOut, cp, negCurveString, i
            end
            if zzNp1 > tr_size^2
                zOut = project_to_boundary_with_coefs(z, d, tr_size,
                                                      zz, zd, dd)
                return zOut, cp, boundaryString, i
            end
            # z = zNp1 # z + alpha*d
            copyto!(z, zNp1)
                
            r += alpha * hess_vec_func(d)
            # Pr = P \ r
            ldiv!(Pr, P, r)
            rPrNp1 = dot(r, Pr)
              
            if dot(r, r) < cg_tol_sq
                return z, cp, interiorString, i
            end
            beta = rPrNp1 / rPr
            rPr = rPrNp1
            # d = -Pr + beta*d
            copyto!(d, -Pr + beta * d)
        
            zz = zzNp1
            zd, dd = cg_inner_products(alpha, beta, zd, dd, rPr, z, d)
        end
    end

    return z, cp, interiorString * "_", settings.max_cg_iters
end

struct TrustRegionSolverGPU{
    V,
    ObjectiveCache,
    Precond,
    Settings,
    Timer
}
    cauchy_point::CauchyPoint{V, Timer}
    dogleg_step::DogLegStep{Timer}
    model_problem::TrustRegionModelProblem{V, Timer}
    objective_cache::ObjectiveCache
    preconditioner::Precond
    settings::Settings
    timer::Timer
    use_warm_start::Bool
    verbose::Bool
end

function TrustRegionSolverGPU(
    objective_cache, p;
    preconditioner = CholeskyPreconditioner,
    settings = TrustRegionSolverSettings(),
    timer = TimerOutput(),
    use_warm_start = false,
    verbose = true
)
    cauchy_point = CauchyPoint(create_unknowns(objective_cache), timer)
    dogleg_step = DogLegStep(timer)
    model_problem = TrustRegionModelProblem(
        create_unknowns(objective_cache),
        create_unknowns(objective_cache),
        create_unknowns(objective_cache),
        create_unknowns(objective_cache),
        timer
    )
    precond = preconditioner(objective_cache, p, timer)
    return TrustRegionSolverGPU(
        cauchy_point, dogleg_step, model_problem,
        objective_cache, precond, 
        settings, timer, use_warm_start, verbose
    )
end

function is_converged(solver::TrustRegionSolverGPU, objective, x, realO, modelO, realRes, modelRes, cgIters, trSize)
    gg = dot(realRes, realRes)
    if gg < solver.settings.tol^2
        modelResNorm = norm(modelRes)
        realResNorm = sqrt(gg)
        print_min_banner(
            solver,
            realO, modelO,
            realResNorm,
            modelResNorm,
            cgIters,
            trSize,
            interiorString,
            true
        )
    
        if solver.verbose
            @info "Converged"
            println("") # a bit of output formatting
        end
    
        if solver.settings.check_stability
            @assert false "hook this up"
            # objective.check_stability(x)
        end    
        return true
    end
    return false
end

function print_banner(solver::TrustRegionSolverGPU)
    if solver.verbose
        @info "=============================================================================================================================="
        @info "Objective        Model Objective  Residual        Model Residual      CG    TR size"
        @info "=============================================================================================================================="
    end
end
  
function print_min_banner(
    solver::TrustRegionSolverGPU,
    objective, modelObjective, res, modelRes, 
    cgIters, trSize, onBoundary, willAccept
)
    if solver.settings.debug_info == false
        return
    end

    if willAccept
        will_accept = true
    else
        will_accept = false
    end

    if solver.verbose
        str = @sprintf "% 1.6e    % 1.6e    %1.6e    %1.6e    %6i    %1.6e    %s    %s" objective modelObjective res modelRes cgIters trSize onBoundary will_accept
        @info str
    end
end

function solve!(solver::TrustRegionSolverGPU, Uu, p)
    @timeit solver.timer "TrustRegionSolver - solve!" begin
        # unpack solver settings
        settings = solver.settings
        tr_size = settings.tr_size

        if solver.use_warm_start
            # TODO implemt GPU compatable warm start
        end

        # calculate initial objective and gradient
        o = value(solver.objective_cache, Uu, p)
        g = gradient(solver.objective_cache, Uu, p)
        g_norm = norm(g)

        if solver.verbose
            @info @sprintf "Initial objective = %1.6e" sum(o)
            @info @sprintf "Initial residual  = %1.6e" g_norm
        end

        # preconditioner
        update_preconditioner!(solver.preconditioner, solver.objective_cache, Uu, p)
        P = solver.preconditioner#.preconditioner

        if is_converged(solver, value, Uu, 0.0, 0.0, g, g, 0, tr_size)
            return nothing
        end

        cumulative_cg_iters = 0
        print_banner(solver)

        # begin loop over trust region iterations
        for i in 1:settings.max_trust_iters
            if i > 1 && i % 25 == 0
                print_banner(solver)
            end

            # setup functions
            if settings.use_incremental_objective
                @assert false "not implemented yet"
            else
                increment_objective = d -> value(solver.objective_cache, Uu + d, p) - o
            end

            hess_vec_func = v -> hvp(solver.objective_cache, Uu, p, v)

            # TODO need to fix below
            # K = hessian!(solver.preconditioner.assembler, solver.objective, x, p)
            # mult_by_approx_hessian = v -> (K + 0.0 * I) * v
            # mult_by_approx_hessian = v -> K \ v
            # mult_by_approx_hessian = v -> hvp(solver.objective, x, p, v)
            mult_by_approx_hessian = v -> v

            # calculate cauchy point
            # cp_norm_sq = calculate!(solver.cauchy_point, g, hess_vec_func, mult_by_approx_hessian, tr_size)
            cp_norm_sq = calculate!(
                solver.cauchy_point, 
                solver.objective_cache, Uu, p, g,
                mult_by_approx_hessian, tr_size
            )

            # minimize model problem
            if cp_norm_sq >= tr_size * tr_size
                @warn "unpreconditioned gradient cauchy point outside trust region at dist = $(sqrt.(cauchyPointNormSquared))"
                # cauchyPoint *= (tr_size / sqrt(cp_norm_sq))
                copyto!(solver.cauchy_point.cauchy_point, (tr_size / sqrt(cp_norm_sq)) * solver.cauchy_point.cauchy_point)
                cp_norm_sq = tr_size * tr_size
                q_newton_point = solver.cauchy_point.cauchy_point
                stepType = boundaryString
                cgIters = 1
            else
                # @assert false "Implement model problem"
                q_newton_point, _, step_type, cg_iters = 
                    minimize!(solver.model_problem, Uu, g, P, hess_vec_func, tr_size, solver.settings)
            end

            cumulative_cg_iters += cg_iters
            tr_size_used = tr_size
            happy_about_tr_size = false

            while !happy_about_tr_size
                d = dogleg_step(solver.cauchy_point.cauchy_point, q_newton_point, tr_size, mult_by_approx_hessian)
                Jd = hess_vec_func(d)
                dJd = dot(d, Jd)
                model_objective = dot(g, d) + 0.5 * dJd
                real_objective = increment_objective(d)
                y = Uu + d
                gy = gradient(solver.objective_cache, y, p)

                if is_converged(
                    solver,
                    solver.objective_cache, y, real_objective, model_objective,
                    gy, g + Jd, cg_iters, tr_size_used
                )
                    # Uu .= y
                    copyto!(Uu, y)
                    return nothing
                end

                model_improve = -model_objective
                real_improve = -real_objective

                rho = real_improve / model_improve

                if model_objective > 0
                    @warn "Found a positive model objective increase.  Debug if you see this."
                    @warn "modelObjective = $model_objective"
                    @warn "realObjective = $real_objective"
                    rho = real_improve / -model_improve
                end

                if !(rho >= settings.η_2)  # write it this way to handle NaNs
                    tr_size *= settings.t_1
                elseif rho > settings.η_3 && is_on_boundary(step_type)
                    tr_size *= settings.t_2
                end

                model_res = g + Jd
                model_res_norm = norm(model_res)
                real_res_norm = norm(gy)

                will_accept = rho >= settings.η_1 || 
                     (rho >= -0 && realResNorm <= g_norm)
                print_min_banner(
                    solver,
                    real_objective, model_objective,
                    real_res_norm, model_res_norm,
                    cg_iters, tr_size_used, step_type,
                    will_accept,
                    # settings
                )

                if will_accept
                    # x = y
                    # g = gy
                    copyto!(Uu, y)
                    copyto!(g, gy)
                    o = value(solver.objective_cache, Uu, p)
                    g_norm = real_res_norm
                    tried_new_precond = false
                    happy_about_tr_size = true
                    # if callback
          
                    # end
                  else
                    # set these for output
                    # trust region will continue to strink until we find a solution on the boundary
                    step_type = boundaryString
                    cg_iters = 0
                end

                if cg_iters >= settings.max_cg_iters || 
                    cumulative_cg_iters >= settings.max_cumulative_cg_iters
                    # objective.update_precond(x)
                    update_preconditioner!(solver.preconditioner, solver.objective_cache, Uu, p; verbose=solver.verbose)
                    P = solver.preconditioner
                    cumulative_cg_iters = 0
                end

                tr_size_used = tr_size

                if tr_size < settings.min_tr_size
        
                    if !tried_new_precond
                        @info "The trust region is too small, updating precond and trying again."
                        update_preconditioner!(solver.preconditioner, solver.objective, Uu, p)
                        P = solver.preconditioner.preconditioner
                        cumulative_cg_iters = 0
                        tried_new_precond = true
                        happy_about_tr_size = true
                        tr_size = settings.tr_size                    
                    else
                        @warn "The trust region is still too small.  Accepting, but be careful."
                        return nothing
                    end
                end
            end
        end
        @assert false "Finish implementing TrustRegionSolverGPU"
    end

    @assert false "reached maximum tr iterations"
end
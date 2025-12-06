# TODO
# make not global ASAP
const negCurveString = "neg curve"
const boundaryString = "boundary"
const interiorString = "interior"


struct CauchyPoint{RV <: AbstractArray{<:Number, 1}}
    cauchy_point::RV
    timer::TimerOutput
end

function calculate!(
    cauchy_point::CauchyPoint, 
    objective_cache,
    Uu,
    p,
    g, 
    hess_vec_func,
    mult_by_approx_hessian,
    tr_size
)
    @timeit cauchy_point.timer "Cauchy Point" begin
        cp = cauchy_point.cauchy_point
        # Kg = hvp(objective_cache, Uu, p, g)
        Kg = hvp(objective_cache, Uu, g, p)
        # Kg = hess_vec_func(g)
        gKg = dot(g, Kg)
        if gKg > 0
            alpha = -dot(g, g) / gKg
            # @time copyto!(cp, alpha .* g)
            cp .= alpha .* g
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

struct DogLegStep{RV <: AbstractArray{<:Number, 1}}
    d::RV
    result::RV
    timer::TimerOutput
end

function preconditioned_project_to_boundary!(r, z, d, tr_size, zz, mult_by_approx_hessian)
    Pd = mult_by_approx_hessian(d)
    dd = dot(d, Pd)
    zd = dot(z, Pd)
    tau = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
    r .= z .+ tau .* d
    return nothing
end

function solve!(dogleg::DogLegStep, cp, newton_p, tr_size, mat_mul)
    @timeit dogleg.timer "DogLegStep - solve!" begin
        cp = cp.cauchy_point
        cc = dot(cp, mat_mul(cp))
        nn = dot(newton_p, mat_mul(newton_p))
        tt = tr_size * tr_size
        
        if cc >= tt #return cauchy point if it extends outside the tr
            dogleg.result .= cp .* sqrt(tt / cc)
            return dogleg.result
        end
    
        if cc > nn # return cauchy point?  seems the preconditioner was not accurate?
            @warn "cp outside newton, preconditioner likely inaccurate"
            return cp
        end
    
        if nn > tt # on the dogleg (we have nn >= cc, and tt >= cc)
            dogleg.d .= newton_p .- cp
            preconditioned_project_to_boundary!(
                dogleg.result,
                cp,
                dogleg.d,
                tr_size,
                cc,
                mat_mul
            )
            return dogleg.result
        end
    end
    copyto!(dogleg.result, newton_p)
    return dogleg.result
end

struct TrustRegionModelProblem{RV <: AbstractArray{<:Number, 1}}
    cp::RV
    d::RV
    Pr::RV
    z::RV
    zNp1::RV
    zOut::RV
    timer::TimerOutput
end

function cg_inner_products_preconditioned(alpha, beta, zd, dd, rPr, z, d)
    # recurrence formulas from Gould et al. doi:10.1137/S1052623497322735
    zd = beta * (zd + alpha * dd)
    dd = rPr + beta * beta * dd
    return zd, dd
end
  
function cg_inner_products_unpreconditioned(alpha, beta, zd, dd, rPr, z, d)
    zd = dot(z, d)
    dd = dot(d, d)
    return zd, dd
end
  
function preconditioned_project_to_boundary(z, d, trSize, zz, mult_by_approx_hessian)
    # find tau s.t. (z + tau*d)^2 = trSize^2
    Pd = mult_by_approx_hessian(d)
    dd = dot(d, Pd)
    zd = dot(z, Pd)
    tau = (sqrt((trSize^2 - zz) * dd + zd^2) - zd) / dd
    # is it ever better to choose the - sqrt() branch?
    return z + tau * d
end
  
function project_to_boundary_with_coefs!(r, z, d, tr_size, zz, zd, dd)
    tau = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
    r .= z .+ tau .* d
    return nothing
end

function update_step_length_squared(alpha, zz, zd, dd)
    return zz + 2 * alpha * zd + alpha * alpha * dd
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
        zNp1 = model_problem.zNp1

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
        copyto!(cp, -Pr)
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
                
            # zNp1 = z + alpha*d
            zNp1 .= z .+ alpha .* d
            zzNp1 = update_step_length_squared(alpha, zz, zd, dd)
        
            if curvature <= 0
                # zOut = project_to_boundary_with_coefs(z, d, tr_size,
                #                                       zz, zd, dd)
                zOut = model_problem.zOut
                project_to_boundary_with_coefs!(zOut, z, d, tr_size, zz, zd, dd)
                return zOut, cp, negCurveString, i
            end
            if zzNp1 > tr_size^2
                zOut = model_problem.zOut
                # zOut = project_to_boundary_with_coefs(z, d, tr_size,
                #                                       zz, zd, dd)
                project_to_boundary_with_coefs!(zOut, z, d, tr_size, zz, zd, dd)
                return zOut, cp, boundaryString, i
            end
            # z = zNp1 # z + alpha*d
            copyto!(z, zNp1)
            
            r += alpha * hess_vec_func(d)
            # r .= r .+ alpha .* hess_vec_func(d)
            
            # Pr = P \ r
            ldiv!(Pr, P, r)
            rPrNp1 = dot(r, Pr)
              
            if dot(r, r) < cg_tol_sq
                return z, cp, interiorString, i
            end
            beta = rPrNp1 / rPr
            rPr = rPrNp1
            # d = -Pr + beta*d
            # copyto!(d, -Pr + beta * d)
            d .= -Pr .+ beta .* d
        
            zz = zzNp1
            zd, dd = cg_inner_products(alpha, beta, zd, dd, rPr, z, d)
        end
    end

    return z, cp, interiorString * "_", settings.max_cg_iters
end

"""
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct TrustRegionSolverSettings #<: NonlinearSolverSettings
    t_1::Float64                                  = 0.25
    t_2::Float64                                  = 1.75
    η_1::Float64                                  = 1e-10
    η_2::Float64                                  = 0.1
    η_3::Float64                                  = 0.5
    max_trust_iters::Int64                        = 100
    tol::Float64                                  = 1.e-8
    max_cg_iters::Int64                           = 50
    max_cumulative_cg_iters::Int64                = 1000
    cg_tol::Float64                               = 0.2e-9
    cg_inexact_solve_ratio::Float64               = 1e-5
    tr_size::Float64                              = 2.0
    min_tr_size::Float64                           = 1e-8
    check_stability::Bool                         = false
    use_preconditioned_inner_product_for_cg::Bool = true
    use_incremental_objective::Bool               = false
    debug_info::Bool                              = true
    over_iter::Int64                              = 0
    verbose::Bool                                 = true
end

function Base.show(io::IO, settings::TrustRegionSolverSettings)
  print(io, "  TrustRegionSolverSettings\n")
  for name in fieldnames(typeof(settings))
    print(io, "        ", rpad(name, 39), " = ", getfield(settings, name), "\n")
  end
  print(io, "\n")
end

function is_on_boundary(stepType)
    return stepType==boundaryString || stepType==negCurveString
end

struct TrustRegionSolver{
    RV <: AbstractArray{<:Number, 1},
    ObjectiveCache,
    Precond,
    WS <: Union{Nothing, WarmStart}
}
    cauchy_point::CauchyPoint{RV}
    dogleg_step::DogLegStep{RV}
    model_problem::TrustRegionModelProblem{RV}
    objective_cache::ObjectiveCache
    preconditioner::Precond
    settings::TrustRegionSolverSettings
    timer::TimerOutput
    warm_start::WS
end

function TrustRegionSolver(
    objective_cache, parameters;
    preconditioner = CholeskyPreconditioner,
    timer = TimerOutput(),
    use_warm_start = false,
    kwargs...
)
    @timeit timer "TrustRegionSolver - setup" begin
        settings = TrustRegionSolverSettings(; kwargs...)
        cauchy_point = CauchyPoint(
            create_unknowns(objective_cache), 
            timer
        )
        dogleg_step = DogLegStep(
            create_unknowns(objective_cache),
            create_unknowns(objective_cache),
            timer
        )
        model_problem = TrustRegionModelProblem(
            create_unknowns(objective_cache),
            create_unknowns(objective_cache),
            create_unknowns(objective_cache),
            create_unknowns(objective_cache),
            create_unknowns(objective_cache),
            create_unknowns(objective_cache),
            timer
        )
        precond = preconditioner(
            objective_cache, parameters, timer
        )

        if use_warm_start
            warm_start = WarmStart(
                objective_cache,
                parameters,
                timer
            )
        else
            warm_start = nothing
        end
    end

    return TrustRegionSolver(
        cauchy_point, dogleg_step, model_problem,
        objective_cache, precond, 
        settings, timer, warm_start
    )
end

function is_converged(solver::TrustRegionSolver, objective, x, realO, modelO, realRes, modelRes, cgIters, trSize)
    @timeit solver.timer "convergence check" begin
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
        
            if solver.settings.verbose
                @info "Converged"
                println("") # a bit of output formatting
            end
        
            if solver.settings.check_stability
                @assert false "hook this up"
                # objective.check_stability(x)
            end    
            return true
        end
    end
    return false
end

function print_banner(solver::TrustRegionSolver)
    if solver.settings.verbose
        @info "=============================================================================================================================="
        @info "Objective        Model Objective  Residual        Model Residual      CG    TR size"
        @info "=============================================================================================================================="
    end
end
  
function print_min_banner(
    solver::TrustRegionSolver,
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

    if solver.settings.verbose
        str = @sprintf "% 1.6e    % 1.6e    %1.6e    %1.6e    %6i    %1.6e    %s    %s" objective modelObjective res modelRes cgIters trSize onBoundary will_accept
        @info str
    end
end

function solve!(solver::TrustRegionSolver, Uu, p)
    @timeit solver.timer "TrustRegionSolver - solve!" begin
        # unpack solver settings
        settings = solver.settings
        tr_size = settings.tr_size


        if solver.warm_start !== nothing
            # preconditioner
            # update_preconditioner!(solver.preconditioner, solver.objective_cache, Uu, p)
            # P = solver.preconditioner#.preconditioner
        
            # TODO implemt GPU compatable warm start
            solve!(
                solver.warm_start, solver.objective_cache, 
                Uu, p; 
                # P=P,
                verbose=solver.settings.verbose
            )
        end

        # calculate initial objective and gradient
        o = value(solver.objective_cache, Uu, p)
        g = gradient(solver.objective_cache, Uu, p)
        g_norm = norm(g)

        if solver.settings.verbose
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

            # hess_vec_func = v -> hvp(solver.objective_cache, Uu, p, v)
            hess_vec_func = v -> hvp(solver.objective_cache, Uu, v, p)

            # TODO need to fix below
            # mult_by_approx_hessian = v -> hvp(solver.objective_cache, Uu, p, v)
            mult_by_approx_hessian = v -> v

            # calculate cauchy point
            cp_norm_sq = calculate!(
                solver.cauchy_point, 
                solver.objective_cache, Uu, p, g,
                hess_vec_func,
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
                # d = dogleg_step(solver.cauchy_point.cauchy_point, q_newton_point, tr_size, mult_by_approx_hessian)
                d = solve!(solver.dogleg_step, solver.cauchy_point, q_newton_point, tr_size, mult_by_approx_hessian)
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
                    update_preconditioner!(solver.preconditioner, solver.objective_cache, Uu, p; verbose=solver.settings.verbose)
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
    end

    @assert false "reached maximum tr iterations"
end

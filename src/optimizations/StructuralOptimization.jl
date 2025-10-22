# TODO lots of room to pre-allocate and re-size objects
# to save on allocations
# this will require a lot of work thought to basically
# put re-init methods all cruical containers

struct StructuralOptimization{
    C, O, Sens, Sim, T,
    # C, O, Sim, T,
    F, DFDPs, DFDX, P, DXDPs
} <: AbstractOptimization{O}
    # geometry::G
    csm_file::C
    objective::O
    sensitivity::Sens
    simulation::Sim
    timer::T
    #
    f::F
    dfdp::DFDPs
    dfdX::DFDX
    p::P
    dXdps::DXDPs
end

function forward_problem! end
function gradient_and_value end
function gradient_check! end
function optimize! end
function update! end
function value end

# function StructuralOptimization(csm_file, objective, sim)
#     sens = ESPSensitivity(csm_file)
#     mesh = UnstructuredMesh(sim.mesh_file) # TODO we can probably get info in another way
#     timer = TimerOutput()

#     # scratch arrays
#     f = zeros(1)
#     dfdp = zeros(length(sens.dxdps))
#     dfdX = similar(mesh.nodal_coords)
#     fill!(dfdX, zero(eltype(dfdX)))

#     ndims = size(mesh.nodal_coords, 1)
#     syms = (:X, :Y)
#     if ndims == 3
#         syms = (syms..., :Z)
#     end

#     # TODO pre-allocate as Vector{H1Field{blah}}
#     dXdps = []
#     for temp in sens.dxdps
#         n_nodes = length(temp) ÷ 3
#         temp_r = reshape(temp, 3, n_nodes)[1:ndims, :]
#         push!(dXdps, H1Field(temp_r, syms))
#     end
#     dXdps = convert(Vector{typeof(dXdps[1])}, dXdps)

#     return StructuralOptimization(
#         csm_file, objective, sens, sim, timer,
#         f,
#         dfdp,
#         dfdX,
#         sens.tesselation.geometry.design_parameter_values,
#         dXdps
#     )
# end

# function forward_problem!(opt::StructuralOptimization, design_params)
#     update!(opt, design_params)
#     sens = opt.sensitivity

#     mesh_file = split(opt.csm_file, ".csm")[1] * "_temp.exo"

#     sim_objective = QuadratureLevelObjective(energy, residual, stiffness)
#     objective_cache = QuadratureLevelObjectiveCache(sim_objective, opt.simulation)
#     Uu = create_unknowns(objective_cache)
#     p = parameters(objective_cache)

#     # TODO need to reset time. Probably a better way to force reset all paramters
#     fill!(p.times.time_current, zero(eltype(p.times.time_current)))

#     assembler = objective_cache.sim_cache.assembler

#     displ = assembler.dof.H1_vars[1]
#     solver = Cthonios.TrustRegionSolver(objective_cache, p, opt.timer; use_warm_start=true, verbose=false)
#     mesh = UnstructuredMesh(mesh_file)
#     pp = PostProcessor(mesh, "output.exo", displ)
#     Cthonios.evolve!(objective_cache.sim_cache, solver, pp, Uu, p)
#     close(pp)

#     # @info "Finished running forward problem"

#     storage_vals = map(similar, values(assembler.scalar_quadarature_storage))
#     storage = NamedTuple{keys(assembler.scalar_quadarature_storage)}(storage_vals)

#     qoi = ScalarQOIExtractor(
#         objective_cache.sim_cache.assembler,
#         energy, sum, storage
#     )

#     qoi_gradient_and_value!(opt.f, opt.dfdX, qoi, Uu, p)
    
#     # update sensitivities
#     for (n, sens) in enumerate(opt.dXdps)
#         opt.dfdp[n] = dot(opt.dfdX.vals, sens.vals)
#     end

#     return sum(opt.f)
# end

# function gradient_and_value(opt::StructuralOptimization, x, g)
#     forward_problem!(opt, x)
#     copyto!(g, opt.dfdp)
#     return value(opt)
# end

# # TODO move to tests
# function gradient_check!(opt::StructuralOptimization, p, step)
#     dfdp = copy(opt.dfdp)
#     d = rand(Float64, length(dfdp))
#     d = d / norm(d)
#     dfdp_fd = (forward_problem!(opt, p .+ step * d) .- forward_problem!(opt, p)) / step
#     return dfdp, dfdp_fd
# end

# function optimize!(opt::StructuralOptimization)
#     x0 = copy(opt.p)
#     p = nothing
#     lb = opt.sensitivity.tesselation.geometry.design_parameters_lbs
#     ub = opt.sensitivity.tesselation.geometry.design_parameters_ubs

#     func = (x, g) -> gradient_and_value(opt, x, g)
    
#     # nl_opt = NLopt.Opt(:LD_MMA, length(x0))
#     nl_opt = NLopt.Opt(:LD_LBFGS, length(x0))
#     NLopt.max_objective!(nl_opt, func)
#     NLopt.lower_bounds!(nl_opt, lb)
#     NLopt.upper_bounds!(nl_opt, ub)
#     NLopt.xtol_rel!(nl_opt, 1e-8)
#     min_f, min_x, ret = NLopt.optimize(nl_opt, x0)
#     # min_f, min_x, ret = NLopt.optimize
# end

# function update!(opt::StructuralOptimization, p_new)
#     # ESP.update_design_parameters!(opt.sensitivity, opt.p)
#     copyto!(opt.p, p_new)
#     update_design_parameters!(opt.sensitivity, p_new)
#     mesh_file = split(opt.csm_file, ".csm")[1] * "_temp.exo"
#     mesh = UnstructuredMesh(mesh_file)
#     ndims = size(mesh.nodal_coords, 1)
#     resize!(opt.dfdX.vals, length(mesh.nodal_coords.vals))
#     for (n, temp) in enumerate(opt.sensitivity.dxdps)
#         n_nodes = length(temp) ÷ 3
#         temp_r = vec(reshape(temp, 3, n_nodes)[1:ndims, :])
#         # temp_r = vec(reshape(temp, n_nodes, 3)[:, 1:ndims])
#         resize!(opt.dXdps[n].vals, length(temp_r))
#         opt.dXdps[n].vals .= temp_r
#     end
#     return nothing
# end

# function value(opt::StructuralOptimization)
#     return sum(opt.f)
# end

module CthoniosESPExt

import Cthonios: ScalarQOIExtractor
import Cthonios: StructuralOptimization
using Cthonios
using EngineeringSketchPadWrapper
using FiniteElementContainers
using LinearAlgebra
using NLopt
using TimerOutputs

function Cthonios.StructuralOptimization(csm_file, objective, sim)
    sens = ESPSensitivity(csm_file)
    mesh = UnstructuredMesh(sim.mesh_file) # TODO we can probably get info in another way
    timer = TimerOutput()

    # scratch arrays
    f = zeros(1)
    dfdp = zeros(length(sens.dxdps))
    dfdX = similar(mesh.nodal_coords)
    fill!(dfdX, zero(eltype(dfdX)))

    ndims = size(mesh.nodal_coords, 1)

    # TODO pre-allocate as Vector{H1Field{blah}}
    dXdps = []
    for temp in sens.dxdps
        n_nodes = length(temp) ÷ 3
        temp_r = reshape(temp, 3, n_nodes)[1:ndims, :]
        push!(dXdps, H1Field(temp_r))
    end
    dXdps = convert(Vector{typeof(dXdps[1])}, dXdps)

    return StructuralOptimization(
        csm_file, objective, sens, sim, timer,
        f,
        dfdp,
        dfdX,
        sens.tesselation.geometry.design_parameter_values,
        dXdps
    )
end

function Cthonios.forward_problem!(opt::StructuralOptimization, design_params)
    @info "New forward problem with parameters $design_params"
    Cthonios.update!(opt, design_params)
    sens = opt.sensitivity

    solver = x -> Cthonios.TrustRegionSolverGPU(x; use_warm_start=true, verbose=false)
    objective_cache = Cthonios.run!(opt.simulation, Cthonios.QuasiStaticObjectiveCache, solver)
    assembler = objective_cache.assembler
    p = objective_cache.parameters
    fill!(p.times.time_current, zero(eltype(p.times.time_current)))

    storage_vals = map(similar, values(objective_cache.assembler.scalar_quadrature_storage))
    storage = NamedTuple{keys(assembler.scalar_quadrature_storage)}(storage_vals)

    qoi = ScalarQOIExtractor(
        objective_cache.assembler,
        energy, sum, storage # TODO make not hardcoded to energy
    )
    Uu = objective_cache.solution
    p = objective_cache.parameters

    Cthonios.qoi_gradient_and_value!(opt.f, opt.dfdX, qoi, Uu, p)
    
    # update sensitivities
    for (n, sens) in enumerate(opt.dXdps)
        opt.dfdp[n] = dot(opt.dfdX.data, sens.data)
    end

    return sum(opt.f)
end

function Cthonios.gradient_and_value(opt::StructuralOptimization, x, g)
    Cthonios.forward_problem!(opt, x)
    copyto!(g, opt.dfdp)
    return Cthonios.value(opt)
end

# TODO move to tests
function Cthonios.gradient_check!(opt::StructuralOptimization, p, step)
    Cthonios.forward_problem!(opt, p)
    dfdp = copy(opt.dfdp)
    d = rand(Float64, length(dfdp))
    d = d / norm(d)
    dfdp_fd = (Cthonios.forward_problem!(opt, p .+ step * d) - Cthonios.forward_problem!(opt, p)) / step
    # dfdp_fd = (Cthonios.forward_problem!(opt, p .+ step) .- Cthonios.forward_problem!(opt, p)) / step

    return dfdp, dfdp_fd
end

function Cthonios.optimize!(opt::StructuralOptimization)
    x0 = copy(opt.p)
    lb = opt.sensitivity.tesselation.geometry.design_parameters_lbs
    ub = opt.sensitivity.tesselation.geometry.design_parameters_ubs

    func = (x, g) -> Cthonios.gradient_and_value(opt, x, g)
    
    # nl_opt = NLopt.Opt(:LD_MMA, length(x0))
    nl_opt = NLopt.Opt(:LD_MMA, length(x0))
    NLopt.min_objective!(nl_opt, func)
    NLopt.lower_bounds!(nl_opt, lb)
    NLopt.upper_bounds!(nl_opt, ub)
    NLopt.xtol_rel!(nl_opt, 1e-8)
    min_f, min_x, ret = NLopt.optimize(nl_opt, x0)
end

function Cthonios.update!(opt::StructuralOptimization, p_new)
    copyto!(opt.p, p_new)
    update_design_parameters!(opt.sensitivity, p_new)
    mesh_file = split(opt.csm_file, ".csm")[1] * "_temp.exo"
    mesh = UnstructuredMesh(mesh_file)
    ndims = size(mesh.nodal_coords, 1)
    resize!(opt.dfdX.data, length(mesh.nodal_coords.data))
    for (n, temp) in enumerate(opt.sensitivity.dxdps)
        n_nodes = length(temp) ÷ 3
        temp_r = vec(reshape(temp, 3, n_nodes)[1:ndims, :])
        resize!(opt.dXdps[n].data, length(temp_r))
        # opt.dXdps[n].data .= temp_r
        copyto!(opt.dXdps[n].data, temp_r)
    end
    return nothing
end

function Cthonios.value(opt::StructuralOptimization)
    return sum(opt.f)
end

end # module

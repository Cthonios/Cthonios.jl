struct PostProcessor{
    ExoPP <: FiniteElementContainers.PostProcessor,
    MO, # TODO type to material output
    MV
}
    exodus_pp::ExoPP
    mat_outputs::MO
    mat_vars::MV
end

function PostProcessor(objective_cache, sim)
    mesh = UnstructuredMesh(sim.mesh_file)
    V_q = FunctionSpace(mesh, L2QuadratureField, Lagrange)
    disp_var = objective_cache.assembler.dof.var

    mat_outputs, mat_vars = create_material_output(objective_cache, V_q, StandardMaterialOutput{Float64})
    update_material_output!(mat_outputs, objective_cache)

    exodus_pp = FiniteElementContainers.PostProcessor(
        mesh, sim.output_file, 
        disp_var, mat_vars...
    )
    pp = PostProcessor(exodus_pp, mat_outputs, mat_vars)
    _post_process_common!(pp, objective_cache, 1)

    return pp
end

function close(pp::PostProcessor)
    FiniteElementContainers.close(pp.exodus_pp)
end

function post_process(
    pp::PostProcessor,
    objective_cache, n
)
    _post_process_common!(pp, objective_cache, n)
end

function _post_process_common!(pp, objective_cache, n)
    mat_outputs, mat_vars = pp.mat_outputs, pp.mat_vars

    U_names = names(assembler(objective_cache).dof.var)
    # U = objective_cache.solution
    p = objective_cache.parameters
    U = p.h1_field
    update_material_output!(mat_outputs, objective_cache)

    write_times(pp.exodus_pp, n, sum(p.times.time_current))
    write_field(pp.exodus_pp, n, U_names, U)
    # write_field(pp, n, )

    for (block, val) in pairs(mat_outputs)
        write_field(pp.exodus_pp, n, String(block), "algorithmic_energy", val.algorithmic_energy)
        write_field(pp.exodus_pp, n, String(block), "cauchy_stress", val.cauchy_stress)
        write_field(pp.exodus_pp, n, String(block), "displacement_gradient", val.displacement_gradient)
    end
end

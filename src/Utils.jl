# TODO does below not play nice with enzyme due to the temp value of field
# being passed to assembler since this is a quadrature storage?
function assemble_scalar_for_ad!(storage, assembler, Uu, p, func)
    # fill!(storage, zero(eltype(storage)))
    fspace = FiniteElementContainers.function_space(assembler.dof)
    t = FiniteElementContainers.current_time(p.times)
    Δt = FiniteElementContainers.time_step(p.times)
    update_field_dirichlet_bcs!(p.h1_field, p.dirichlet_bcs)
    update_field_unknowns!(p.h1_field, assembler.dof, Uu)
    for (b, (field, conns, block_physics, state_old, state_new, props)) in enumerate(zip(
        values(storage),
        values(fspace.elem_conns), 
        values(p.physics),
        values(p.state_old), values(p.state_new),
        values(p.properties)
    ))
        ref_fe = values(fspace.ref_fes)[b]
        backend = FiniteElementContainers._check_backends(assembler, p.h1_field, p.h1_coords, state_old, state_new, conns)
        # FiniteElementContainers._assemble_block_scalar!(
        FiniteElementContainers._assemble_block_quadrature_quantity!(
            field, block_physics, ref_fe, 
            p.h1_field, p.h1_field_old, p.h1_coords, state_old, state_new, props, t, Δt,
            # conns, b, residual,
            # conns, b, func,
            conns, b, energy,
            backend
        )
        # display(field)
    end
end

# function assemble_vector_for_ad!(storage, assembler, Uu, p, func)
#     fill!(storage, zero(eltype(storage)))
#     fspace = FiniteElementContainers.function_space(assembler.dof)
#     t = FiniteElementContainers.current_time(p.times)
#     Δt = FiniteElementContainers.time_step(p.times)
#     update_field_dirichlet_bcs!(p.h1_field, p.dirichlet_bcs)
#     update_field_unknowns!(p.h1_field, assembler.dof, Uu)
#     for (b, (conns, block_physics, state_old, state_new, props)) in enumerate(zip(
#         values(fspace.elem_conns), 
#         values(p.physics),
#         values(p.state_old), values(p.state_new),
#         values(p.properties)
#     ))
#         ref_fe = values(fspace.ref_fes)[b]
#         # backend = FiniteElementContainers._check_backends(assembler, p.h1_field, p.h1_coords, state_old, state_new, conns)
#         backend = KA.get_backend(p.h1_field)
#         FiniteElementContainers._assemble_block_vector!(
#             storage, block_physics, ref_fe, 
#             p.h1_field, p.h1_field_old, p.h1_coords, state_old, state_new, props, t, Δt,
#             # conns, b, residual,
#             conns, b, func,
#             backend
#         )
#     end
# end

# function create_properties(physics, props)
#     @assert length(physics) == length(props)
#     for (k1, k2) in zip(keys(physics), keys(props))
#         @assert k1 == k2
#     end
#     props = map(
#         (x, y) -> FiniteElementContainers.create_properties(x, y), 
#         values(physics), values(props)
#     )
#     props = NamedTuple{keys(physics)}(props)
#     return props
# end

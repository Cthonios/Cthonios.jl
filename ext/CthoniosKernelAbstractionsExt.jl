module CthoniosKernelAbstractionsExt

import KernelAbstractions as KA
using Cthonios
using FiniteElementContainers
using KernelAbstractions

# Nonlinear solver stuff

@kernel function update_bcs_kernel!(U, coords, t, bc)
  I = @index(Global)
  X = @views coords[:, bc.nodes[I]]
  U[bc.dof, I] = bc.func(X, t)
  nothing
end

function Cthonios.update_bcs!(U, coords, t, bcs, backend::Cthonios.KernelAbstractionsBackend{<:CPU})
  kernel = update_bcs_kernel!(backend.backend)
  for bc in bcs
    kernel(U, coords, t, bc, ndrange=length(bc.nodes))
  end
end

@kernel function residual_kernel!(R, section, U, X)
  ND, NN       = num_dimensions(section), num_nodes_per_element(section)
  model, props = section.model, section.props
  state        = SVector{0, Float64}()

  I = @index(Global)
  conn = connectivity(section, I)
  dof_conn = dof_connectivity(section, I)
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views vec(X[:, conn]))
  U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views vec(U[:, conn]))
  # R_el = zero(SMatrix{ND, NN, Float64, ND * NN})
  R_el = zero(SVector{ND * NN, eltype(U)}) # Unitful issue here

  for q in 1:num_q_points(section)
    ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
    J    = X_el * ∇N_ξ
    J_inv = inv(J)
    ∇N_X = (J_inv * ∇N_ξ')'
    JxW  = det(J) * ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe, q)
    ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
    F_q  = ∇u_q + one(∇u_q)
    P_q  = ConstitutiveModels.pk1_stress(model, props, F_q, state)
    P_v  = FiniteElementContainers.extract_stress(section.formulation, P_q) 
    G    = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)
    R_el   = R_el + JxW * G * P_v
  end
  FiniteElementContainers.assemble_residual!(R, R_el, dof_conn)
end

function Cthonios.residual!(R, section::Cthonios.Section, U::NodalField, X::NodalField, backend::Cthonios.KernelAbstractionsBackend{CPU})
  kernel = residual_kernel!(backend.backend)
  kernel(R, section, U, X)
  nothing
end

# function Cthonios.solve!(
#   solver::Cthonios.NewtonSolver, domain::Cthonios.QuasiStaticDomain,
#   ka_backend::Cthonios.KernelAbstractionsBackend{CPU}
# )

#   # TODO probably need to have an adapt for the solver prior to calling solve
#   # TODO just assume that's the case for now.
#   Uu, ΔUu = solver.Uu, solver.ΔUu
#   U = FiniteElementContainers.create_fields(ka_backend.backend, domain.dof)
#   R = FiniteElementContainers.create_fields(ka_backend.backend, domain.dof)
#   K = domain.assembler.K

#   # TODO need an actual kernel for this in FEMContainers
#   Cthonios.update_bcs!(U, domain.coords, domain.time.current_time, domain.bcs, ka_backend)

#   # TODO need an update_stiffness kernel

#   @assert false
# end
# @kernel function energy_kernel()

# end

# function Cthonios.energy()

# end

end # module
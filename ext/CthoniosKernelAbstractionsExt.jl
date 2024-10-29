module CthoniosKernelAbstractionsExt

using Atomix
using ComponentArrays
using Cthonios
using FiniteElementContainers
using KernelAbstractions
using StaticArrays

# utilities
@kernel function zero_grad_kernel!(vals)
  I = @index(Global)
  vals[I] = zero(eltype(vals))
end

function zero_grad!(vals)
  backend = get_backend(vals.vals)
  kernel = zero_grad_kernel!(backend)
  kernel(vals, ndrange=length(vals.vals))
  return nothing
end

# domain kernels
@kernel function update_dirichlet_vals!(Ubc, @Const(bc), @Const(X), @Const(t), @Const(prev_total))
  I = @index(Global)
  index = prev_total + I
  node = bc.nodes[I]
  X_temp = @views X[:, node]
  val = bc.func(X_temp, t)
  Ubc[index] = val
end

# currently inefficient but works
# TODO this will break if you change up the bc locations
function Cthonios.update_dirichlet_vals!(Ubc, domain::Domain, X, t, backend::Backend)
  kernel = update_dirichlet_vals!(backend)
  # t = t.current_time
  t = Cthonios.current_time(t)
  prev_total = 0
  for bc in domain.dirichlet_bcs
    kernel(Ubc, bc, X, t, prev_total, ndrange=length(bc.nodes))
    prev_total = prev_total + length(bc.nodes)
  end
  return nothing
end

@kernel function update_field_bcs!(U, @Const(dirichlet_dofs), @Const(Ubc))
  I = @index(Global)
  U[dirichlet_dofs[I]] = Ubc[I]
end

function Cthonios.update_field_bcs!(U, domain::Domain, Ubc, backend::Backend)
  kernel = update_field_bcs!(backend)
  kernel(U, domain.dirichlet_dofs, Ubc, ndrange=length(Ubc))
end

# move this kernel to FiniteElementContainers KA ext
@kernel function update_field_unknowns!(U, @Const(dof), @Const(Uu))
  I = @index(Global)
  U[dof.unknown_dofs[I]] = Uu[I]
end 

# move this kernel to FiniteElementContainers KA ext
function Cthonios.update_field_unknowns!(U, domain::Domain, Uu, backend::Backend)
  kernel = update_field_unknowns!(backend)
  kernel(U, domain.dof, Uu, ndrange=length(Uu))
  return nothing
end

# objective methods
@kernel function gradient!(vals, @Const(section), @Const(U), @Const(X), @Const(props))
  # Q, E = @index(Global, NTuple)
  E = @index(Global)

  ND, NN, NP, NS = size(section)
  NF = Cthonios.num_fields(section)
  fspace, physics = section.fspace, section.physics

  conn = connectivity(fspace, E)
  dof_conn = dof_connectivity(fspace, E)
  U_el = Cthonios.element_fields(section, U, dof_conn)
  f_el = zeros(SVector{NF * NN, eltype(vals)})
  # TODO note, this will break for any problem where
  # NF != ND
  for Q in 1:FiniteElementContainers.num_q_points(fspace)
    X_el = Cthonios.element_fields(section, X, dof_conn)
    # TODO clean up below a lot
    N    = FiniteElementContainers.shape_function_values(fspace, Q)
    ∇N_ξ = FiniteElementContainers.shape_function_gradients(fspace, Q)
    ∇N_X = FiniteElementContainers.map_shape_function_gradients(X_el, ∇N_ξ)
    JxW  = FiniteElementContainers.volume(X_el, ∇N_ξ) * FiniteElementContainers.quadrature_weights(fspace, Q)
    X_q  = X_el * N
    interps = FiniteElementContainers.Interpolants(X_q, N, ∇N_X, JxW)

    f_el += Cthonios.gradient(physics, interps, U_el, props)
  end

  # assemble
  for i in axes(dof_conn, 1)
    Atomix.@atomic vals.vals[dof_conn[i]] += f_el[i]
  end
end

function Cthonios.gradient!(vals, o::Objective, Uu, p, backend::Backend)
  zero_grad!(vals)
  Cthonios.update_field_bcs!(p.U, o.domain, p.Ubc, backend)
  Cthonios.update_field_unknowns!(p.U, o.domain, Uu, backend)
  kernel = gradient!(backend)
  for (name, sec) in pairs(o.domain.sections)
    NE = FiniteElementContainers.num_elements(sec.fspace)
    sec_props = @view p.props[name]
    kernel(vals, sec, p.U, p.X, sec_props, ndrange=(NE,))
  end
  return nothing
end

function Cthonios.gradient(o::Objective, Uu, p, backend::Backend)
  vals = FiniteElementContainers.create_fields(o.domain.dof, backend)
  Cthonios.gradient!(vals, o, Uu, p, backend)
  return view(vals, o.domain.dof.unknown_dofs)
end

# function Cthonios.hessian(o::Objective, Uu, p, backend::Backend)

# end

# TODO currently assumes fixed properties per section
@kernel function objective!(vals, @Const(section), @Const(U), @Const(X), @Const(props))
  Q, E = @index(Global, NTuple)

  ND, NN, NP, NS = size(section)
  NF = Cthonios.num_fields(section)
  fspace, physics = section.fspace, section.physics

  conn = connectivity(fspace, E)
  dof_conn = dof_connectivity(fspace, E)
  U_el = Cthonios.element_fields(section, U, dof_conn)
  # TODO note, this will break for any problem where
  # NF != ND
  X_el = Cthonios.element_fields(section, X, dof_conn)
  # TODO clean up below a lot
  N    = FiniteElementContainers.shape_function_values(fspace, Q)
  ∇N_ξ = FiniteElementContainers.shape_function_gradients(fspace, Q)
  ∇N_X = FiniteElementContainers.map_shape_function_gradients(X_el, ∇N_ξ)
  JxW  = FiniteElementContainers.volume(X_el, ∇N_ξ) * FiniteElementContainers.quadrature_weights(fspace, Q)
  X_q  = X_el * N
  interps = FiniteElementContainers.Interpolants(X_q, N, ∇N_X, JxW)

  vals[Q, E] = Cthonios.energy(physics, interps, U_el, props)
end

# remove this
@kernel function objective_v2!(vals, @Const(section), @Const(U), @Const(X), @Const(props))
  # Q, E = @index(Global, NTuple)
  E = @index(Global)

  ND, NN, NP, NS = size(section)
  NF = Cthonios.num_fields(section)
  fspace, physics = section.fspace, section.physics

  conn = connectivity(fspace, E)
  dof_conn = dof_connectivity(fspace, E)
  U_el = Cthonios.element_fields(section, U, dof_conn)
  # TODO note, this will break for any problem where
  # NF != ND
  for Q in 1:FiniteElementContainers.num_q_points(fspace)
    X_el = Cthonios.element_fields(section, X, dof_conn)
    # TODO clean up below a lot
    N    = FiniteElementContainers.shape_function_values(fspace, Q)
    ∇N_ξ = FiniteElementContainers.shape_function_gradients(fspace, Q)
    ∇N_X = FiniteElementContainers.map_shape_function_gradients(X_el, ∇N_ξ)
    JxW  = FiniteElementContainers.volume(X_el, ∇N_ξ) * FiniteElementContainers.quadrature_weights(fspace, Q)
    X_q  = X_el * N
    interps = FiniteElementContainers.Interpolants(X_q, N, ∇N_X, JxW)

    vals[Q, E] = Cthonios.energy(physics, interps, U_el, props)
  end
end

function Cthonios.objective!(o::Objective, Uu, p, backend::Backend)
  Cthonios.update_field_bcs!(p.U, o.domain, p.Ubc, backend)
  Cthonios.update_field_unknowns!(p.U, o.domain, Uu, backend)
  kernel = objective!(backend)
  for (name, sec) in pairs(o.domain.sections)
    NQ = FiniteElementContainers.num_q_points(sec.fspace)
    NE = FiniteElementContainers.num_elements(sec.fspace)
    sec_vals = @view p.q_vals_scratch[name]
    sec_props = @view p.props[name]
    kernel(sec_vals, sec, p.U, p.X, sec_props, ndrange=(NQ, NE))
    # kernel(sec_vals, sec, p.U, p.X, sec_props, ndrange=(NE,))
  end
  # return sum(p.q_vals_scratch)
  return nothing
end

function Cthonios.objective(o::Objective, Uu, p, backend)
  Cthonios.objective!(o, Uu, p, backend)
  return sum(p.q_vals_scratch)
end

end # module

using Cthonios
# using Enzyme
using ForwardDiff
using FiniteElementContainers
using LinearAlgebra
using SparseArrays
using StaticArrays

# file management
mesh_file = Base.source_dir() * "/mesh.g"

# global setup
times = ConstantTimeStepper(0.0, 1.0, 0.0125)
n_dofs = 1

# functions
zero_bc(x, t) = 0.0
one_bc(x, t) = 1.0
func(x) = 2π^2 * sin(2π * x[1]) * sin(2π * x[2])

# bcs
dirichlet_bcs = DirichletBC[
  DirichletBC("unnamed_nset_1", [1], zero_bc)
  DirichletBC("nset_left", [1], zero_bc)
  DirichletBC("unnamed_nset_3", [1], zero_bc)
  DirichletBC("nset_right", [1], zero_bc)
]
neumann_bcs = [
]

# sections
sections = [
  TotalLagrangeSection("unnamed_block_1", Poisson(func), 2)
]

# domain setup
domain = Domain(
  mesh_file, times, 
  sections, dirichlet_bcs, neumann_bcs, 
  n_dofs
)
@time Cthonios.update_unknown_dofs!(domain)
@time Cthonios.update_bcs!(domain)

objective = Objective(
  domain,
  Cthonios.energy!,
  Cthonios.gradient!,
  Cthonios.hessian!
)
# # solver = DirectSolver(domain)
# solver = Cthonios.NewtonSolver(DirectSolver, domain)
solver = Cthonios.NewtonSolver(DirectSolver, objective)

# # # testing stuff out
Uu = Cthonios.create_unknowns(objective.domain)
# # residual = Cthonios.create_fields(domain)
# # hess = solver.assembler
p = nothing
# # en = [0.0]

@time Cthonios.update_fields!(objective.domain, Uu)

Cthonios.solve!(solver, objective, Uu, p)
# Uu

# @time Cthonios.solve!(solver, objective, Uu, p)
# @time Cthonios.solve!(solver, objective, Uu, p)


# @time Cthonios.energy!(en, objective, Uu, p)
# @time Cthonios.energy!(en, objective, Uu, p)

# @time Cthonios.gradient!(residual, objective, Uu, p)
# @time Cthonios.gradient!(residual, objective, Uu, p)

# @time Cthonios.hessian!(hess, objective, Uu, p)
# @time Cthonios.hessian!(hess, domain, Uu, p)
# SparseArrays.sparse!(hess) |> Symmetric

# # NN = size(domain.X, 1)
# section = domain.sections[1]
# fspace = domain.sections[1].fspace
# conn = connectivity(fspace, 1)
# dof_conn = dof_connectivity(fspace, 1)
# X = Cthonios.element_coordinates(section, domain.X, conn)
# NN = size(X, 2)
# # U = element_fields(section, domain.U, dof_conn)
# U = SVector{NN * 1, eltype(domain.U)}(@views domain.U[dof_conn])
# state = SVector{0, Float64}()
# props = SVector{0, Float64}()
# ref_fe = fspace.ref_fe
# N      = FiniteElementContainers.ReferenceFiniteElements.shape_function_values(ref_fe, 1)
# ∇N_ξ   = FiniteElementContainers.ReferenceFiniteElements.shape_function_gradients(ref_fe, 1)
# J      = (X * ∇N_ξ)'
# J_inv  = inv(J)
# ∇N_X   = (J_inv * ∇N_ξ')'
# JxW    = det(J) * FiniteElementContainers.ReferenceFiniteElements.quadrature_weights(ref_fe, 1)


# @time Cthonios.element_iterator!(en, Cthonios.quadrature_energy, domain, Uu, p)
# @time Cthonios.element_iterator!(en, Cthonios.quadrature_energy, domain, Uu, p)

# # @time Cthonios.quadrature_energy(section.physics, U, state, X, 0.0, props, N, ∇N_X, JxW, ∇N_X)
# # @time Cthonios.quadrature_energy(section.physics, U, state, X, 0.0, props, N, ∇N_X, JxW, ∇N_X)

# @time Cthonios.quadrature_gradient(section.physics, U, state, X, 0.0, props, N, ∇N_X, JxW, ∇N_X)
# @time Cthonios.quadrature_gradient(section.physics, U, state, X, 0.0, props, N, ∇N_X, JxW, ∇N_X)

# @time Cthonios.quadrature_hessian(section.physics, U, state, X, 0.0, props, N, ∇N_X, JxW, ∇N_X)
# @time Cthonios.quadrature_hessian(section.physics, U, state, X, 0.0, props, N, ∇N_X, JxW, ∇N_X)

# # @time ForwardDiff.gradient(
# #   z -> Cthonios.quadrature_energy(section.physics, z, state, X, 0.0, props, N, ∇N_X, JxW, ∇N_X), 
# #   U
# # )
# # @time ForwardDiff.hessian(
# #   z -> Cthonios.quadrature_energy(section.physics, z, state, X, 0.0, props, N, ∇N_X, JxW, ∇N_X), 
# #   U
# # )

# # @time autodiff(
# #   Reverse, Cthonios.quadrature_energy, Active,
# #   Const(section.physics), Active(U), Const(state), Const(X), 
# #   Const(0.0), Const(props), Const(N), Const(∇N_X), 
# #   Const(JxW), Const(∇N_X)
# # )
# # @time autodiff(
# #   Reverse, Cthonios.quadrature_energy, Active,
# #   Const(section), Active(U), Const(state), Active(X), 
# #   Const(0.0), Const(props), Const(N), Const(∇N_X), 
# #   Const(JxW), Const(∇N_X)
# # )

# # Enzyme.hessian()

# # en
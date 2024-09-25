module CthoniosKernelAbstractionsExt

using ComponentArrays
using Cthonios
using FiniteElementContainers
using KernelAbstractions

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
function Cthonios.update_dirichlet_vals!(Ubc, domain::Domain, X, t, backend)
  kernel = update_dirichlet_vals!(backend)
  t = t.current_time
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

function Cthonios.update_field_bcs!(U, domain::Domain, Ubc, backend)
  kernel = update_field_bcs!(backend)
  kernel(U, domain.dirichlet_dofs, Ubc, ndrange=length(Ubc))
end

@kernel function update_field_unknowns!(U, @Const(dof), @Const(Uu))
  I = @index(Global)
  U[dof.unknown_dofs[I]] = Uu[I]
end 

function Cthonios.update_field_unknowns!(U, domain::Domain, Uu, backend)
  kernel = update_field_unknowns!(backend)
  kernel(U, domain.dof, Uu, ndrange=length(Uu))
  return nothing
end

# objective methods
@kernel function objective!(vals, @Const(section), @Const(Uu), p, sec_name)
  Q, E = @index(Global, NTuple)
  vals[sec_name][Q, E] = 1.0
end

function Cthonios.objective(o::Objective, Uu, p, backend)
  kernel = objective!(backend)
  for (name, sec) in pairs(o.domain.sections)
    NQ = FiniteElementContainers.num_q_points(sec.fspace)
    NE = FiniteElementContainers.num_elements(sec.fspace)
    kernel(p.q_vals_scratch, sec, Uu, p, name, ndrange=(NQ, NE))
  end
  return sum(vals)
end

end # module

module CthoniosKernelAbstractionsExt

using Cthonios
using KernelAbstractions

@kernel function update_field_bcs!(U, @Const(dirichlet_dofs), @Const(Ubc))
  I = @index(Global)
  U[dirichlet_dofs[I]] = Ubc[I]
end

function Cthonios.update_field_bcs!(U, domain::Domain, Ubc, backend)
  kernel = update_field_bcs!(backend)
  kernel(U, domain.dirichlet_dofs, Ubc, ndrange=length(Ubc))
end

@kernel function update_field_unknowns!(U, @Const(domain), @Const(Uu))
  I = @index(Global)
  U[domain.dof.unknown_dofs[I]] = Uu[I]
end 

function Cthonios.update_field_unknowns!(U, domain::Domain, Uu, backend)
  kernel = update_field_unknowns!(backend)
  kernel(U, domain, Uu, ndrange=length(Uu))
  return nothing
end



end # module

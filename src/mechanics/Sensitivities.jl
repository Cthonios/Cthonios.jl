# Enzyme.API.runtimeActivity!(true)
function Cthonios.energy_gradient(mode::EnzymeCore.ReverseMode, solver, domain::Cthonios.QuasiStaticDomain)
  backend = CPU()

  # preallocate tangent cache
  # tangent_cache = similar(domain.domain_cache)

  # unpack stuff
  Uu = solver.Uu
  X, U, state, props, Π, Πs, V = unpack(domain.domain_cache)
  # dX, dU, dstate, dprops, dΠ, dΠs, dV = unpack(tangent_cache)
  sections = domain.sections
  assembler = solver.linear_solver.assembler
  # # section = sections.section_1

  # # tangent seeding
  # dX .= 0.0
  # dU .= 0.0
  # dstate .= 0.0
  # dprops .= 0.0
  # dΠs .= 0.0
  # dΠ .= 1.0

  # update fields here
  update_fields!(U, domain, X, Uu)

  # now do sensitivies with respect to the whole gambit, even bc nodes
  # energy!(Π, Πs, sections, U, state, props, X, backend)
  # # @time internal_force!()

  # force_kernel! = internal_force_kernel!(backend)
  # kernel_iterator!(assembler, force_kernel!, sections, U, state, props, X)

  stiffness_kernel_temp! = internal_force_and_stiffness_kernel!(backend)
  kernel_iterator!(assembler, stiffness_kernel_temp!, sections, U, state, props, X)

  # @time autodiff(
  #   Reverse, energy!,
  #   Duplicated(Π, dΠ),
  #   Duplicated(Πs, dΠs),
  #   Const(sections),
  #   Duplicated(U, dU),
  #   Duplicated(state, dstate),
  #   Duplicated(props, dprops),
  #   Duplicated(X, dX),
  #   Const(backend)
  # )

  # dX
end
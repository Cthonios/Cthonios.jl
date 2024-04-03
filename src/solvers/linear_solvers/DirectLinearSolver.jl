# CURRENTLY NOT HOOKED UP TO ANYTHING
"""
$(TYPEDFIELDS)
"""
struct DirectLinearSolverSettings{F1, F2} <: AbstractLinearSolverSettings
  factorization_method::F1
  factorization_method!::F2
end

"""
$(TYPEDSIGNATURES)
"""
function DirectLinearSolverSettings(input_settings::Dict{Symbol, Any})
  method_name = input_settings[Symbol("factorization method")]
  in_place_method_name = method_name * "_factorize!"
  factorization_method = eval(Meta.parse(method_name))
  factorization_method! = eval(Meta.parse(in_place_method_name))
  return DirectLinearSolverSettings(factorization_method, factorization_method!)
end

"""
$(TYPEDFIELDS)
"""
struct DirectLinearSolver{Settings, Assembler, Factorization} <: AbstractLinearSolver{Settings, Assembler}
  settings::Settings
  assembler::Assembler
  factorization::Factorization
end

function Base.show(io::IO, solver::DirectLinearSolver)
  print(io, "DirectLinearSolver\n")
  # print(io, "  $(solver.settings)\n")
  print(io, "      $(typeof(solver.assembler).name.name)\n")
  print(io, "      $(typeof(solver.factorization))\n")
end

# TODO maybe add some settings
"""
$(TYPEDSIGNATURES)
"""
function DirectLinearSolver(input_settings::Dict{Symbol, Any}, domain::QuasiStaticDomain, backend)
  settings = DirectLinearSolverSettings(input_settings)
  assembler = StaticAssembler(domain.dof, map(x -> x.fspace, values(domain.sections)))

  # need to update dofs prior to factorization
  update_unknown_dofs!(assembler, domain)

  # also need to assemble once
  Δt = domain.domain_cache.time.Δt 
  U = domain.domain_cache.U
  # Uu = FiniteElementContainers.create_unknowns(domain.dof) # TODO remove this allocation
  Uu = domain.domain_cache.Uu
  X = domain.domain_cache.X
  update_fields!(U, domain, X, Uu)
  
  sections = domain.sections
  state_old = domain.domain_cache.state_old
  state_new = domain.domain_cache.state_new
  props = domain.domain_cache.props
  stiffness!(assembler, state_new, sections, Δt, X, U, props, state_old, backend)

  # setup matrix to setup a factorization
  # TODO eventually we might have matrix free stuff
  # TODO so this is probably not general
  K = SparseArrays.sparse!(assembler)
  factorization = ldl(K)

  return DirectLinearSolver(settings, assembler, factorization)
end

"""
$(TYPEDSIGNATURES)
"""
function update_unknown_dofs!(solver::DirectLinearSolver, domain::QuasiStaticDomain)
  update_unknown_dofs!(solver.assembler, domain)
end

"""
$(TYPEDSIGNATURES)
"""
function solve!(Uu, solver::DirectLinearSolver, domain::QuasiStaticDomain, common::CthoniosCommon)
  K = SparseArrays.sparse!(solver.assembler)
  # R = @views domain.domain_cache.f[domain.dof.unknown_dofs]
  R = @views solver.assembler.residuals[domain.dof.unknown_dofs]
  @timeit timer(common) "Factorization" begin
    # solver.settings.factorization_method!(K, solver.factorization)
    ldl_factorize!(K, solver.factorization)
  end
  @timeit timer(common) "Back substitution" begin
    ldiv!(Uu, solver.factorization, R)
  end
  return nothing
end

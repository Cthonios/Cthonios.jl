struct DirectLinearSolverSettings{F1, F2} <: AbstractLinearSolverSettings
  factorization_method::F1
  factorization_method!::F2
end

function DirectLinearSolverSettings(input_settings::Dict{Symbol, Any})
  method_name = input_settings[Symbol("factorization method")]
  in_place_method_name = method_name * "!"
  factorization_method = eval(Meta.parse(method_name))
  factorization_method! = eval(Meta.parse(in_place_method_name))
  return DirectLinearSolverSettings(factorization_method, factorization_method!)
end

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
function DirectLinearSolver(input_settings::Dict{Symbol, Any}, domain::QuasiStaticDomain)
  settings = DirectLinearSolverSettings(input_settings)
  assembler = StaticAssembler(domain.dof, map(x -> x.fspace, values(domain.sections)))

  # need to update dofs prior to factorization
  update_unknown_dofs!(assembler, domain)

  # also need to assemble once
  U = domain.domain_cache.U
  Uu = FiniteElementContainers.create_unknowns(domain.dof) # TODO remove this allocation
  X = domain.coords
  update_fields!(U, domain, X, Uu)
  
  sections = domain.sections
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  stiffness!(assembler, U, sections, state, props, X)

  # setup matrix to setup a factorization
  # TODO eventually we might have matrix free stuff
  # TODO so this is probably not general
  K = SparseArrays.sparse!(assembler) |> Symmetric
  factorization = settings.factorization_method(K)

  return DirectLinearSolver(settings, assembler, factorization)
end

function solve!(Uu, solver::DirectLinearSolver, domain::QuasiStaticDomain, common::CthoniosCommon)
  K = SparseArrays.sparse!(solver.assembler) |> Symmetric
  R = solver.assembler.residuals[domain.dof.unknown_dofs]
  @timeit timer(common) "Factorization" begin
    solver.settings.factorization_method!(solver.factorization, K)
  end
  @timeit timer(common) "Back substitution" begin
    Uu .= solver.factorization \ R
  end
  return nothing
end

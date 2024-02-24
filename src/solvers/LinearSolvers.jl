abstract type AbstractLinearSolverSettings end
abstract type AbstractLinearSolver{Settings, Assembler} end

function solve! end

function update_residual!(
  solver::AbstractLinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}

  state = domain.domain_cache.state
  props = domain.domain_cache.props
  internal_force!(solver.assembler.residuals, U, domain, Uu, state, props, domain.coords)
  return nothing
end

function update_stiffness!(
  solver::AbstractLinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}

  state = domain.domain_cache.state
  props = domain.domain_cache.props
  stiffness!(solver.assembler, U, domain, Uu, state, props, domain.coords)
  return nothing
end

function update_residual_and_stiffness!(
  solver::AbstractLinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}

  state = domain.domain_cache.state
  props = domain.domain_cache.props
  internal_force_and_stiffness!(solver.assembler, U, domain, Uu, state, props, domain.coords)
  return nothing
end

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
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  stiffness!(assembler, U, domain, Uu, state, props, domain.coords)

  # setup matrix to setup a factorization
  # TODO eventually we might have matrix free stuff
  # TODO so this is probably not general
  K = SparseArrays.sparse!(assembler) |> Symmetric
  factorization = settings.factorization_method(K)

  return DirectLinearSolver(settings, assembler, factorization)
end

function update_unknown_dofs!(assembler, d::QuasiStaticDomain)
  # update the dofs
  FiniteElementContainers.update_unknown_dofs!(d.dof, d.bc_dofs)
  FiniteElementContainers.update_unknown_dofs!(assembler, d.dof, map(x -> x.fspace, d.sections), d.bc_dofs)
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

function setup_linear_solver(input_settings, domain)
  type = Meta.parse(input_settings[:type])
  return eval(type)(input_settings, domain)
end

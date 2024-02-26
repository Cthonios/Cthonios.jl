abstract type AbstractLinearSolverSettings end
abstract type AbstractLinearSolver{Settings, Assembler} end

function solve! end

function update_unknown_dofs!(assembler::StaticAssembler, d::QuasiStaticDomain)
  # update the dofs
  FiniteElementContainers.update_unknown_dofs!(d.dof, d.bc_dofs)
  FiniteElementContainers.update_unknown_dofs!(assembler, d.dof, map(x -> x.fspace, d.sections), d.bc_dofs)
end

# function update_unknown_dofs!(assembler::MatrixFreeStaticAssembler, d::QuasiStaticDomain)
#   # update the dofs
#   FiniteElementContainers.update_unknown_dofs!(d.dof, d.bc_dofs)
#   # FiniteElementContainers.update_unknown_dofs!(assembler, d.dof, map(x -> x.fspace, d.sections), d.bc_dofs)
# end

# DirectLinearSolver
include("DirectLinearSolver.jl")
include("IterativeLinearSolver.jl")

# general setup
function setup_linear_solver(input_settings, domain)
  type = Meta.parse(input_settings[:type])
  return eval(type)(input_settings, domain)
end
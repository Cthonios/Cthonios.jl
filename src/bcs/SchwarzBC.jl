struct SchwarzDirchletBC{D, S, B} <: AbstractBCInput
  dirichlet_bc::D
  coupled_subsim::S
  coupled_subsim_block_name::B
end
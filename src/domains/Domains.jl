abstract type AbstractDomainCache end

abstract type AbstractDomain{
  Dof,
  Funcs,
  BCNodes,
  BCDofs,
  BCFuncIDs,
  Sections
} end

include("QuasiStaticDomain.jl")

function setup_domain(input_settings)
  new_section("Domains")
  
end
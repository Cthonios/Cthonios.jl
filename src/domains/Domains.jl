abstract type AbstractDomainCache end

abstract type AbstractDomain{
  Dof,
  Funcs,
  BCNodes,
  BCDofs,
  BCFuncIDs,
  Sections,
  Time,
  DomainCache
} end

include("QuasiStaticDomain.jl")

function setup_domain(input_settings)
  new_section("Domains")
  
end
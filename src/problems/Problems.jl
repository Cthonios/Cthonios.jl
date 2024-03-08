function get_domain_input_settings(inputs::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any} 
  return inputs[:domain]
end

function get_mesh_file_name(inputs::D)::String where D <: Dict{Symbol, Any}
  return inputs[:domain][:mesh][Symbol("file name")]
end

function get_output_file_name(inputs::D)::String where D <: Dict{Symbol, Any}
  return inputs[:results][Symbol("output file name")]
end

function get_output_nodal_fields(inputs::D)::Vector{String} where D <: Dict{Symbol, Any}
  return inputs[:results][Symbol("nodal fields")]
end

function get_output_element_fields(inputs::D)::Vector{String} where D <: Dict{Symbol, Any}
  if Symbol("element fields") in keys(inputs[:results])
    return inputs[:results][Symbol("element fields")]
  else
    return Vector{String}(undef, 0)
  end
end

function get_output_quadrature_fields(inputs::D)::Vector{String} where D <: Dict{Symbol, Any}
  if Symbol("quadrature fields") in keys(inputs[:results])
    return inputs[:results][Symbol("quadrature fields")]
  else
    return Vector{String}(undef, 0)
  end
end

function get_solver_input_settings(inputs::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any}
  return inputs[:solver]
end

function get_max_properties(domain)
  return maximum(map(x -> ConstitutiveModels.num_properties(x.model), domain.sections))
end

function get_max_state_variables(domain)
  return maximum(map(x -> ConstitutiveModels.num_state_vars(x.model), domain.sections))
end

function get_max_q_points(domain)
  return maximum(map(x -> FiniteElementContainers.num_q_points(x), domain.sections))
end

include("ForwardProblems.jl")
include("SensitivityProblems.jl")

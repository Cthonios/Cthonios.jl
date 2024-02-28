"""
This method is a barrier method for type stability
"""
function load_input_file(input_file::String)::Dict{Symbol, Any}
  YAML.load_file(input_file; dicttype=Dict{Symbol, Any})
end

"""
Barrier for type stability. Returns input_settings[:boundary_conditions]
but also replaces function names with function input settings.
"""
function parse_boundary_conditions(input_settings::D, domain_key::Symbol) where D <: Dict{Symbol, Any}
  @warn "Only supporting displacement boundary conditions at the moment"
  bc_settings = Dict{Symbol, Any}()
  bc_settings[:displacement] = Vector{Dict{Symbol, Any}}(undef, 0)

  for (n, bc) in enumerate(input_settings[:domains][domain_key][Symbol("boundary conditions")][Symbol("displacement")])
    bc[:function] = input_settings[:functions][Symbol(bc[:function])]
    bc[:function][:func_id] = n
    push!(bc_settings[:displacement], bc)
  end
  return bc_settings
end

"""
Barrier for type stability. Simply returns input_settings[:functions]
"""
function parse_functions(input_settings::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any}
  return input_settings[:functions]
end

"""
Barrier for type stability. Simply returns input_settings[Symbol("linear solvers")]
"""
function parse_linear_solvers(input_settings::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any}
  return input_settings[Symbol("linear solvers")]
end

"""
Barrier for type stability. Simply returns input_settings[:materials]
"""
function parse_materials(input_settings::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any}
  return input_settings[:materials]
end

"""
Barrier for type stability. Simply returns input_settings[Symbol("nonlinear solvers")]
"""
function parse_nonlinear_solvers(input_settings::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any}
  return input_settings[Symbol("nonlinear solvers")]
end

"""
Barrier for type stability. Simply returns input_settings[:mesh]
"""
function parse_mesh(input_settings::D, domain_key::Symbol)::Dict{Symbol, String} where D <: Dict{Symbol, Any}
  return input_settings[:domains][domain_key][:mesh]
end 

"""
Barrier for type stability. Returns input_settings[:sections] but
replaces material names with material definitions
"""
function parse_sections(input_settings::D, domain_key::Symbol)::Vector{Dict{Symbol, Any}} where D <: Dict{Symbol, Any}
  section_settings = input_settings[:domains][domain_key][:sections]
  for section in section_settings
    section[:material] = input_settings[:materials][Symbol(section[:material])]
  end

  return section_settings
end

"""
Barrier for type stability. Simple returns input_settings[Symbol("time stepper")]
"""
function parse_time_stepper(input_settings::D, domain_key::Symbol)::Dict{Symbol, Any} where D <: Dict{Symbol, Any}
  return input_settings[:domains][domain_key][Symbol("time stepper")]
end


#############################

"""
Helper to parse one domain
"""
function parse_domains(input_settings::D, key) where D <: Dict{Symbol, Any}
  domain_settings = Dict{Symbol, Any}()
  domain_settings[:mesh] = parse_mesh(input_settings, key)
  domain_settings[:boundary_conditions] = parse_boundary_conditions(input_settings, key)
  domain_settings[:sections] = parse_sections(input_settings, key)
  domain_settings[:time_stepper] = parse_time_stepper(input_settings, key)
  return domain_settings
end

"""
Helper to parse a domain and all settings need to set it up
"""
function parse_domains(input_settings::D) where D <: Dict{Symbol, Any}
  domain_settings = Dict{Symbol, Any}()
  # for key in keys(input_settings[:domains])
  #   domain_settings[key] = Dict{Symbol, Any}()
  #   domain_settings[key][:mesh] = parse_mesh(input_settings, key)
  #   domain_settings[key][:boundary_conditions] = parse_boundary_conditions(input_settings, key)
  #   domain_settings[key][:sections] = parse_sections(input_settings, key)
  #   domain_settings[key][:time_stepper] = parse_time_stepper(input_settings, key)
  # end
  for key in keys(input_settings[:domains])
    domain_settings[key] = parse_domains(input_settings, key)
  end
  return domain_settings
end 

"""
Helper to parse problems
"""
function parse_problems(
  input_settings::D, domain_settings::D, 
  linear_solver_settings::D, nonlinear_solver_settings::D
) where D <: Dict{Symbol, Any}
  problem_settings = Dict{Symbol, Any}()
  for key in keys(input_settings[:problems])
    problem_settings[key] = Dict{Symbol, Any}()
    problem_settings[key][:type] = input_settings[:problems][key][:type]
    problem_settings[key][:domain] = domain_settings[Symbol(input_settings[:problems][key][:domain])]
    problem_settings[key][:solver] = nonlinear_solver_settings[Symbol(input_settings[:problems][key][:solver])]
    problem_settings[key][:results] = input_settings[:problems][key][:results]

    # now swap out the nonlinear solver linear solver settings
    problem_settings[key][:solver][Symbol("linear solver")] =
    linear_solver_settings[Symbol(problem_settings[key][:solver][Symbol("linear solver")])]
  end
  return problem_settings
end

"""
Main entry point for parsing
"""
function parse_input_file(input_file::String)
  # load file
  loaded_settings = load_input_file(input_file)

  # Pre-allocate dict to return
  input_settings = Dict{Symbol, Any}()

  # parse domains scope
  # input_settings[:domains] = parse_domains(loaded_settings)
  domain_settings = parse_domains(loaded_settings)
  linear_solver_settings = parse_linear_solvers(loaded_settings)
  nonlinear_solver_settings = parse_nonlinear_solvers(loaded_settings)

  problem_settings = parse_problems(
    loaded_settings, domain_settings, linear_solver_settings, nonlinear_solver_settings
  )

  input_settings[:domains] = domain_settings
  input_settings[:problems] = problem_settings

  return input_settings
end

function peak_mesh_dims(input_settings::D, key)::Tuple{Int64, Int64} where D <: Dict{Symbol, Any}
  mesh = read_mesh(input_settings[:domains][key])
  return FiniteElementContainers.num_dimensions(mesh) |> Int64,
         FiniteElementContainers.num_nodes(mesh) |> Int64
end


module CthoniosCLIExt

using ArgParse
using ConstitutiveModels
using Cthonios
using FiniteElementContainers
using RuntimeGeneratedFunctions
using YAML

RuntimeGeneratedFunctions.init(@__MODULE__)

function Cthonios.cthonios_main()::Cint
  cli_args = _parse_command_line()
  @assert cli_args["backend"] in ["cpu", "cuda", "rocm"]
  @assert cli_args["input-file"] !== nothing
  @info "cthonios_main"
  @info "Backend    = $(cli_args["backend"])"
  @info "Input file = $(cli_args["input-file"])"
  input_settings = _parse_input_file(cli_args)

  sim_settings = input_settings[:cthonios]
  type = sim_settings[Symbol("simulation type")]

  if type == "SingleDomainSimulation"
    _parse_single_domain_simulation(sim_settings)
  else
    @assert false "unsupported simulation type = $type"
  end

  return Int32(0)
end

function _parse_command_line()::Dict{String, Any}
  settings = ArgParseSettings()
  @add_arg_table! settings begin
    "--backend", "-b"
      arg_type = String
      default = "cpu"
      help = "Backend to use, e.g. cpu, cuda, rocm"
    "--input-file", "-i"
      arg_type = String
      help = "Path to an input file"
      required = true
  end
  return parse_args(settings)
end

function _parse_bcs(settings, type)
  bc_settings = settings[Symbol("$type boundary conditions")]
  # bcs = Dict{Symbol, Any}()
  bcs = []
  for bc in bc_settings
    func = @RuntimeGeneratedFunction(Meta.parse(bc[:function]))
    for field in bc[:fields]
      for sset in bc[:sidesets]
        push!(bcs, DirichletBC(field, sset, func))
      end
    end
  end
  return bcs
end

function _parse_function_spaces(settings, mesh)
  fspace_settings = settings[Symbol("function spaces")]
  fspaces = Dict{Symbol, FunctionSpace}()
  for fspace in fspace_settings
    field_type = eval(Symbol(fspace[Symbol("field type")]))
    interpolation_type = eval(Symbol(fspace[Symbol("interpolation type")]))
    fspace_name = Symbol(fspace[:name])
    fspaces[fspace_name] = FunctionSpace(mesh, field_type, interpolation_type)
  end
  return fspaces
end 

function _parse_input_file(settings)
  input_file = settings["input-file"]
  return YAML.load_file(input_file; dicttype=Dict{Symbol, Any})
end

function _parse_materials(settings)
  material_settings = settings[:materials]
  material_models = Dict{Symbol, Any}()
  material_model_properties = Dict{Symbol, Any}()
  for material in material_settings
    name = Symbol(material[:name])
    material_models[name] = eval(Symbol(material[Symbol("material model")]))
    material_model_properties[name] = material[Symbol("material model properties")]
  end
  return material_models, material_model_properties
end

function _parse_meshes(settings)
  mesh_settings = settings[:mesh]
  type = eval(Symbol(mesh_settings[:type]))
  file_name = mesh_settings[Symbol("file name")]
  return type(file_name)
end

function _parse_nonlinear_solvers(settings)
  solver_settings = settings[Symbol("nonlinear solvers")]
  solvers = Dict{Symbol, Any}()
  solver_settings = Dict{Symbol, Any}()
  for solver in solver_settings
    name = solver_settings[:name]
    type = solver_settings[:type]
    settings = solver_settings[:settings]
    solvers[name] = type
    solver_settings[name] = settings
  end
  return solvers, solver_settings
end

function _parse_objective(settings)
  objective_settings = settings[:objective]
  type = eval(Symbol(objective_settings[:type]))
  value = eval(Symbol(objective_settings[:value]))
  gradient = eval(Symbol(objective_settings[:gradient]))
  hessian = eval(Symbol(objective_settings[:hessian]))
  return type(value, gradient, hessian)
end

function _parse_physics(settings, material_models, material_model_properties)
  physics_settings = settings[:physics]
  type = eval(Symbol(physics_settings[:type]))
  formulation = eval(Symbol(physics_settings[:formulation]))()
  materials = physics_settings[Symbol("material assignment")]
  physics = Dict{Symbol, Any}()
  props = Dict{Symbol, Any}()

  for (material, blocks) in materials
    temp_props = Dict{String, Any}()
    for (name, prop) in material_model_properties[material]
      temp_props[String(name)] = prop
    end

    for block in blocks
      physics[Symbol(block)] = type(formulation, material_models[material]())
      props[Symbol(block)] = temp_props
    end
  end
  physics = NamedTuple(physics)
  props = NamedTuple(props)
  props = Cthonios.create_properties(physics, props)
  return physics, props
end

function _parse_single_domain_simulation(sim_settings)
  mesh = _parse_meshes(sim_settings)
  @info "Mesh:"
  display(mesh)

  fspaces = _parse_function_spaces(sim_settings, mesh)
  @info "Function Spaces:"
  for (name, fspace) in fspaces
    @info "Function Space: $name"
    display(fspace)
  end

  solution_fields = _parse_solution_fields(sim_settings, fspaces)
  @info "Solution Fields:"
  for (name, field) in solution_fields
    @info "Solution Field: $name"
    display(field)
  end

  material_models, material_model_properties = _parse_materials(sim_settings)
  @info "Materials:"
  for (name, model, props) in zip(keys(material_models), values(material_models), values(material_model_properties))
    @info "Material: $name"
    @info "  Model: $model"
    @info "  Properties:"
    for (prop_name, prop) in props
      @info "    $prop_name = $prop"
    end
  end

  physics, props = _parse_physics(sim_settings, material_models, material_model_properties)
  @info "Physics:"
  for (name, physics, props) in zip(keys(physics), values(physics), values(props))
    @info "  Block: $name"
    @info "  Physics: $physics"
    @info "  Properties: $props"
  end

  objective = _parse_objective(sim_settings)
  @info "Objective: $objective"

  # TODO time integrator

  # TODO initial conditions
  dirichlet_bcs = _parse_bcs(sim_settings, "dirichlet")
  @info "Dirichlet Boundary Conditions:"
  for bc in dirichlet_bcs
    @info bc
  end
  dirichlet_bcs = convert(Vector{DirichletBC}, dirichlet_bcs)

  neumann_bcs = _parse_bcs(sim_settings, "neumann")
  @info "Neumann Boundary Conditions:"
  for bc in neumann_bcs
    @info bc
  end
  neumann_bcs = convert(Vector{NeumannBC}, neumann_bcs)


end

function _parse_solution_fields(settings, fspaces)
  field_settings = settings[Symbol("solution fields")]
  fields = Dict{Symbol, Any}()
  for field in field_settings
    field_type = eval(Symbol(field[:type]))
    fspace = fspaces[Symbol(field[Symbol("function space")])]
    field_name = Symbol(field[:name])
    fields[field_name] = field_type(fspace, field_name)
  end
  return fields
end

end # module

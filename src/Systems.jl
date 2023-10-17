include("Variables.jl")


function setup_mesh(input_settings)
  mesh_file = input_settings["mesh"]["file name"]

  str = section_header("Mesh")
  str = str * "Mesh file = $mesh_file \n\n"

  block_names = String[]
  fspace_inputs = input_settings["function spaces"]
  for fspace in fspace_inputs
    push!(block_names, fspace["block"])
  end
  str = str * "Reading blocks:\n"
  for block in block_names
    str = str * "  $block \n"
  end
  str = str * "\n"
  str = str * "Opening exodus database\n"
  exo = ExodusDatabase(mesh_file, "r")
  str = str * "  Reading coordinates\n"
  coordinates = read_coordinates(exo)
  str = str * "  Reading blocks\n"
  # below only needed for efficient implementation
  _, I, B, _ = Exodus.int_and_float_modes(exo.exo)
  blocks = Vector{Block{I, B}}(undef, length(block_names))
  for n in axes(blocks, 1)
    str = str * "    Reading block $(block_names[n])\n"
    blocks[n] = read_set(exo, Block, block_names[n])
  end
  str = str * "\n"
  close(exo)

  str = str * "Mesh successfully read.\n"
  @info str

  return coordinates, blocks
end

function setup_function_spaces(input_settings, coordinates, blocks)
  fspaces_settings = input_settings["function spaces"]
  fspaces = Vector{FunctionSpace}(undef, length(blocks))
  for (n, fspace_settings) in enumerate(fspaces_settings)
    re_name   = fspace_settings["reference finite element"]
    q_degree  = fspace_settings["quadrature degree"]
    el_symbol = Symbol(re_name)
    el_type   = @eval $el_symbol($q_degree)
    re        = ReferenceFE(el_type)
    fspace = FunctionSpace(coordinates, blocks[n], re)
    fspaces[n] = fspace
  end
  return fspaces
end

function setup_variables(input_settings)
  variable_names = input_settings["variables"]
  variable_list  = VariableList(variable_names)
  @show variable_list
end

function System(input_settings)
  coordinates, blocks = setup_mesh(input_settings)
  # variable_list       = V
  setup_variables(input_settings)
  fspaces             = setup_function_spaces(input_settings, coordinates, blocks)
  # dof_manager         = DofManager()
end

"""
Command line interface parser
"""
function cli_parser()
  settings = ArgParseSettings()

  @add_arg_table! settings begin
    "--input-file", "-i"  
      help = "A path to an input file. Currently supported formats are YAML"
      arg_type = String
    "--num-processors", "--np"
      help = "The number of processors to to run the simulation with."
      arg_type = Int
      default = 1
    "--backend"
      help = "The backend to use. Options currently supported are cpu"
      arg_type = String
      default = "cpu"
  end

  return parse_args(settings)
end


# Exceptions for input deck parsing
abstract type CthoniosInputFileParserException <: CthoniosException end

struct MeshInputBlockNotFound <: CthoniosInputFileParserException
end

Base.show(io::IO, ::MeshInputBlockNotFound) = 
print(io, "\nInput file has no mesh block!\n",
      "The syntax for a mesh block is:\n\n",
      "mesh:\n",
      "  file name: <string>\n")

function mesh_input_block_not_found()
  e = MeshInputBlockNotFound()
  @error e
  throw(e)
end 

struct MeshFileNotFound <: CthoniosInputFileParserException
  file_name::String
end

Base.show(io::IO, e::MeshFileNotFound) = print(io, "Mesh file $(e.file_name) not found!\n")

function mesh_file_not_found(file_name::String)
  e = MeshFileNotFound(file_name)
  @error e
  throw(e)
end

struct FunctionSpacesInputBlockNotFound <: CthoniosInputFileParserException
end

Base.show(io::IO, ::FunctionSpacesInputBlockNotFound) = 
print(io, "\nInput file has no function spaces input block!\n",
      "The syntax for a function spaces input block is:\n\n",
      "function spaces:\n",
      "- name:                     <string>\n",
      "  block:                    <string>\n",
      "  quadarture degree:        <int>\n",
      "  reference finite element: <Type<ReferenceFE>>\n",
      "...\n")

function function_spaces_input_block_not_found()
  e = FunctionSpacesInputBlockNotFound()
  @error e
  throw(e)
end

function input_file_parser(input_file::String)
  input_settings = YAML.load_file(input_file)

  # input blocks checking
  if !("mesh" in keys(input_settings))
    mesh_input_block_not_found()
  end

  if !("function spaces" in keys(input_settings))
    function_spaces_input_block_not_found()
  end

  if !isfile(input_settings["mesh"]["file name"])
    mesh_file_not_found(input_settings["mesh"]["file name"])
  end

  # TODO much more to do here
  return input_settings
end


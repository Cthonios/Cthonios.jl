

struct PostProcessor{Out <: ExodusDatabase, V <: NodalField}
  # out_file_name::
  out_file::Out
  scratch_U::V
end

# function PostProcessor(f::FileMesh, out_file::String, dof::DofManager, dims::Int)
"""
Constructor for post processor

closes exodus db by default, needs to be re-opend downstream
this is mainly for debuggin purposes
"""
function PostProcessor(mesh_file::String, out_file::String, dof::DofManager, dims::Int)
  f = FileMesh(ExodusDatabase, mesh_file)
  copy_mesh(f.file_name, out_file)
  Exodus.close(f.mesh_obj)
  # out = FileMesh(ExodusDatabase, out_file)
  out = ExodusDatabase(out_file, "rw")
  if dims == 2
    Exodus.write_names(out, NodalVariable, ["displ_x", "displ_y"])
  elseif dims == 3
    Exodus.write_names(out, NodalVariable, ["displ_x", "displ_y", "displ_z"])
  else
    @assert false "only dim 2, 3 is supported right now"
  end
  # Exodus.close(out)
  U = FiniteElementContainers.create_fields(dof)
  return PostProcessor(out, U)
end

function Base.show(io::IO, pp::PostProcessor)
  print(io, "        PostProcessor\n",
        "          Output file = $(pp.out_file.file_name)\n")
  # printing requested nodal variables
  nodal_vars   = String[]
  element_vars = String[]
  global_vars  = String[]
  for var in keys(pp.out_file.nodal_var_name_dict)
    push!(nodal_vars, var)
  end 
  for var in keys(pp.out_file.element_var_name_dict)
    push!(element_vars, var)
  end
  sort!(nodal_vars)
  print(io, "          Nodal variables\n")
  if length(nodal_vars) > 1
    for var in nodal_vars
      print(io, "            $var\n")
    end
  else
    print(io, "none\n")
  end
  print(io, "          Element variables\n")
  if length(element_vars) > 1
    for var in element_vars
      print(io, "            $var\n")
    end
  else
    print(io, "            none\n")
  end
  print(io, "          Global varialbes\n")
  if length(global_vars) > 1
    for var in global_vars
      print(io, "            $var\n")
    end
  else
    print(io, "            none\n")
  end
end

write_time(p::PostProcessor, n::Int, t::T) where T <: Number = 
Exodus.write_time(p.out_file, n, t)

write_values(p::PostProcessor, type::Type{NodalVariable}, n::Int, name::String, u) = 
Exodus.write_values(p.out_file, type, n, name, u)

close(p::PostProcessor) = Exodus.close(p.out_file)

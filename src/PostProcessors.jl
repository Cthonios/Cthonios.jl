struct PostProcessor{Out <: ExodusDatabase}
  out_file::Out
end

NODAL_FIELDS = Dict{String, String}(
  "displacement"   => "displ",
  "internal force" => "internal_force",
  "dcoordinates"   => "dcoordinates"
)
ELEMENT_FIELDS = Dict{String, String}(
  "properties"  => "properties",
  "dproperties" => "dproperties"
)

function vector_names(base_name::String, n_dims::Int)
  if n_dims == 2
    names = [base_name * "_x", base_name * "_y"]
  elseif n_dims == 3
    names = [base_name * "_x", base_name * "_y", base_name * "_z"]
  end
  return names
end

function variable_length_field_names(base_name::String, n_vars::Int)
  names = String[]
  for n in 1:n_vars
    push!(names, base_name * "_$(lpad(n, 3, '0'))")
  end
  return names
end

"""
Constructor for post processor

closes exodus db by default, needs to be re-opend downstream
this is mainly for debuggin purposes
"""
function PostProcessor(
  # mesh_file::String, 
  domain,
  out_file::String, 
  output_nodal_fields::Vector{String},
  output_element_fields::Vector{String},
  output_quadrature_fields::Vector{String},
  # dims::Int,
  # n_properties::Int, n_state_vars::Int, max_q_points::Int;
  force::Bool = false
)
  # f = FileMesh(ExodusDatabase, mesh_file)
  f = domain.static.mesh
  dims = FiniteElementContainers.num_dimensions(domain.static.mesh) |> Int64

  if isfile(out_file)
    if force
      rm(out_file; force=force)
    end
  end
  copy_mesh(f.file_name, out_file)
  Exodus.close(f.mesh_obj)
  out = ExodusDatabase(out_file, "rw")

  nodal_fields = Vector{String}(undef, 0)
  element_fields = Vector{String}(undef, 0)
  quadrature_fields = Vector{String}(undef, 0)

  # Nodal fields
  for nodal_field in output_nodal_fields
    try
      append!(nodal_fields, vector_names(NODAL_FIELDS[nodal_field], dims))
    catch e
      throw(e)
    end
  end

  # Element fields
  for element_field in output_element_fields
    if element_field == "properties"
      append!(element_fields, variable_length_field_names("properties", n_properties))
    elseif element_field == "dproperties"
      append!(element_fields, variable_length_field_names("dproperties", n_properties))
    else
      @assert false "Unsupported element field requested for output $element_field"
    end
  end

  # Quadrature fields which in Exodus are really element fields
  for quadrature_field in output_quadrature_fields
    if quadrature_field == "state variables"
      for n in 1:n_state_vars
        for q in 1:max_q_points
          push!(element_fields, "state_variables_$(lpad(n, 3, '0'))_$(lpad(q, 3, '0')))")
        end
      end
    else
      @assert false "Unsupport quadrature field requested for output"
    end
  end

  Exodus.write_names(out, NodalVariable, nodal_fields)

  if length(element_fields) > 0
    Exodus.write_names(out, ElementVariable, element_fields)
  end

  # Exodus.close(out)
  return PostProcessor(out)
end

function Base.show(io::IO, pp::PostProcessor)
  print(io, "    PostProcessor\n",
        "      Output file = $(pp.out_file.file_name)\n")
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
  print(io, "      Nodal variables\n")
  if length(nodal_vars) > 1
    for var in nodal_vars
      print(io, "        $var\n")
    end
  else
    print(io, "none\n")
  end
  print(io, "      Element variables\n")
  if length(element_vars) > 1
    for var in element_vars
      print(io, "        $var\n")
    end
  else
    print(io, "        none\n")
  end
  print(io, "      Global varialbes\n")
  if length(global_vars) > 1
    for var in global_vars
      print(io, "        $var\n")
    end
  else
    print(io, "        none\n")
  end
end

write_time(p::PostProcessor, n::Int, t::T) where T <: Number = 
Exodus.write_time(p.out_file, n, t)

write_values(p::PostProcessor, type::Type{NodalVariable}, n::Int, name::String, u) = 
Exodus.write_values(p.out_file, type, n, name, u)

write_values(p::PostProcessor, type::Type{ElementVariable}, n::Int, block_id::Int, name::String, u) = 
Exodus.write_values(p.out_file, type, n, block_id, name, u)

close(p::PostProcessor) = Exodus.close(p.out_file)

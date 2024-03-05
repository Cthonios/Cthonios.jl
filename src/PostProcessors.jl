struct PostProcessor{Out <: ExodusDatabase}
  out_file::Out
end

"""
Constructor for post processor

closes exodus db by default, needs to be re-opend downstream
this is mainly for debuggin purposes
"""
function PostProcessor(
  mesh_file::String, 
  out_file::String, 
  output_nodal_fields::Vector{String},
  output_element_fields::Vector{String},
  output_quadrature_fields::Vector{String},
  dims::Int,
  n_properties::Int, n_state_vars::Int, max_q_points::Int
)
  f = FileMesh(ExodusDatabase, mesh_file)
  copy_mesh(f.file_name, out_file)
  Exodus.close(f.mesh_obj)
  out = ExodusDatabase(out_file, "rw")

  nodal_fields = Vector{String}(undef, 0)
  element_fields = Vector{String}(undef, 0)
  quadrature_fields = Vector{String}(undef, 0)

  # Nodal fields
  for nodal_field in output_nodal_fields
    if nodal_field == "displacement"
      if dims == 2
        append!(nodal_fields, ["displ_x", "displ_y"])
      elseif dims == 3
        append!(nodal_fields, ["displ_x", "displ_y", "displ_z"])
      else
        @assert false "only dim 2, 3 is supported right now"
      end
    elseif nodal_field == "internal force"
      if dims == 2
        append!(nodal_fields, ["internal_force_x", "internal_force_y"])
      elseif dims == 3
        append!(nodal_fields, ["internal_force_x", "internal_force_y", "internal_force_z"])
      else
        @assert false "only dim 2, 3 is supported right now"
      end
    else
      @assert false "Unsupported nodal field $nodal_field"
    end
  end

  # Element fields
  for element_field in output_element_fields
    if element_field == "properties"
      for n in 1:n_properties
        push!(element_fields, "properties_$(lpad(n, 3, '0'))")
      end
    else
      @assert false "Unsupported element field requested for output"
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

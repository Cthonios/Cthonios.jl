

struct PostProcessor{Out <: ExodusDatabase}
  out_file::Out
end

function PostProcessor(f::FileMesh, out_file::String, dims::Int)
  copy_mesh(f.file_name, out_file)
  # out = FileMesh(ExodusDatabase, out_file)
  out = ExodusDatabase(out_file, "rw")
  if dims == 2
    Exodus.write_names(out, NodalVariable, ["displ_x", "displ_y"])
  else
    @assert false "only dim 2 is supported right now"
  end
  return PostProcessor(out)
end

write_time(p::PostProcessor, n::Int, t::T) where T <: Number = 
Exodus.write_time(p.out_file, n, t)

write_values(p::PostProcessor, type::Type{NodalVariable}, n::Int, name::String, u) = 
Exodus.write_values(p.out_file, type, n, name, u)

close(p::PostProcessor) = Exodus.close(p.out_file)

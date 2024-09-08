abstract type AbstractPostProcessor end

struct ExodusPostProcessor{
  Exo <: ExodusDatabase
} <: AbstractPostProcessor
  exo::Exo
  solution_fields::Vector{String}
end

function ExodusPostProcessor(
  mesh_file::String, 
  output_file::String,
  solution_fields::Vector{String}
)
  copy_mesh(mesh_file, output_file)
  exo = ExodusDatabase(output_file, "rw")
  write_names(exo, NodalVariable, solution_fields)
  return ExodusPostProcessor(exo, solution_fields)
end

close(pp::ExodusPostProcessor) = Exodus.close(pp.exo)

function write_fields(pp::ExodusPostProcessor, U, n)
  @assert size(U, 1) == length(pp.solution_fields) "Have size(U, 1) = $(size(U, 1)) and length(pp.solution_fields) = $(length(pp.solution_fields))"

  for (i, name) in enumerate(pp.solution_fields)
    write_values(pp.exo, NodalVariable, n, name, U[i, :])
  end
end

function write_time(pp::ExodusPostProcessor, n, t)
  Exodus.write_time(pp.exo, n, t)
end

export ExodusPostProcessor

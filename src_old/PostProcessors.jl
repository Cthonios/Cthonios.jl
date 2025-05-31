"""
$(TYPEDEF)
Abstract base type for post processors.
"""
abstract type AbstractPostProcessor end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Exodus post processor for writing quantities such
as nodal, element, and quadrature fields to exodus
files via ```Exodus.jl```.
```exo``` - Handle to an open exodus database.
```solution_fields``` - A vector of names for solution fields on the nodes
"""
struct ExodusPostProcessor{
  Exo <: ExodusDatabase
} <: AbstractPostProcessor
  exo::Exo
  solution_fields::Vector{String}
end

"""
$(TYPEDSIGNATURES)
Constructor for ```ExodusPostProcessor```.
```mesh_file``` - A name of a mesh file to copy.
```output_file``` - A name of file to copy the mesh to.
```solution_fields``` - Names for the nodal solution fields.
"""
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

function ExodusPostProcessor(inputs::Dict{Symbol, Any}, mesh_file)
  output_file = inputs[Symbol("output file")]
  field_vars = inputs[Symbol("field variables")]
  return ExodusPostProcessor(mesh_file, output_file, field_vars)
end

"""
$(TYPEDSIGNATURES)
Closes an ```ExodusPostProcessor```
```pp``` - Handle to an open ```ExodusPostProcessor```
"""
close(pp::ExodusPostProcessor) = Exodus.close(pp.exo)

"""
$(TYPEDSIGNATURES)
Writes nodal solution fields stored in ```U``` to 
```pp.exo``` at time step ```n```.
```pp``` - Handle to open ```ExodusPostProcessor```.
```U``` - Nodal solution field stored as a ```NodalField```.
```n``` - Time step.
"""
function write_fields(pp::ExodusPostProcessor, U, n)
  @assert size(U, 1) == length(pp.solution_fields) "Have size(U, 1) = $(size(U, 1)) and length(pp.solution_fields) = $(length(pp.solution_fields))"

  for (i, name) in enumerate(pp.solution_fields)
    write_values(pp.exo, NodalVariable, n, name, U[i, :])
  end
end

"""
$(TYPEDSIGNATURES)
Writes a new time step to ```pp.exo```
```pp.exo``` at time step ```n```.
```pp``` - Handle to open ```ExodusPostProcessor```.
```n``` - Time step.
"""
function write_time(pp::ExodusPostProcessor, n, t)
  Exodus.write_time(pp.exo, n, t)
end

function write_output(pp::ExodusPostProcessor, n, objective, Uu, p)
  @timeit timer(objective) "ExodusPostProcessor - write_output" begin
    update_field_unknowns!(current_solution(p), objective.domain, Uu)
    write_time(pp, n, current_time(p.t))
    write_fields(pp, current_solution(p), n)
  end
end

export ExodusPostProcessor

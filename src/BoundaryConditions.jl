abstract type AbstractBC{Locations, Func} end


# typed so that the return type has consistent type with input
# coordinates. This allows for say units of meters for displacement
# and coordinates and units of seconds for time
struct DisplacementBC{
  D, # dimension of problem
  T1    <: Number, # type for spatial quantity
  T2    <: Number, # type for time quantity
  # Nodes <: AbstractArray{<:Integer, 1},
  Nodes <: NodalField, 
  Func  <: FunctionWrapper{T1, Tuple{SVector{D, T1}, T2}}
} <: AbstractBC{Nodes, Func}

  nodes::Nodes
  dof::Int
  func::Func
end

# TODO make zero of type DistType
function DisplacementBC{D, DistType, TimeType}(
  nodes::Nodes, dof::Int, 
  func::F = (x, t) -> 0.
) where {D, DistType <: Number, TimeType <: Number, 
         Nodes, F <: Function}

  func_wrapper = FunctionWrapper{DistType, Tuple{SVector{D, DistType}, TimeType}}(func)
  return DisplacementBC{D, DistType, TimeType, Nodes, typeof(func_wrapper)}(nodes, dof, func_wrapper)
end

function DisplacementBC{D, DistType, TimeType}(
  nset::N, dof::Int, 
  func::F = (x, t) -> 0.
) where {D, DistType <: Number, TimeType <: Number, N <: NodeSet, F <: Function}

  return DisplacementBC{D, DistType, TimeType}(nset.nodes, dof, func)
end

function DisplacementBC{D, DistType, TimeType}(
  nodes::N, dof::Int,
  func_string::String
) where {D , DistType <: Number, TimeType <: Number, N <: NodalField}

  func_temp = eval(Meta.parse(func_string))
  return DisplacementBC{D, DistType, TimeType}(nodes, dof, func_temp)
end

# # struct EssentialBC{V <: AbstractArray{<:Integer}, F}
# #   nodes::V
# #   dof::Int
# #   func::F
# # end

# struct EssentialBC{
#   D, # dimension of problem
#   T <: 
#   V <: NodalField,
#   F <: FunctionWrapper{, Tuple{<:AbstractArray{<:Number, 1}, <:Number}}
# }
#   nodes::V
#   dof::Int
#   func::F
# end

# # TODO currently only supports floats
# function EssentialBC(nodes::N, dof::Int, func::Function = (x, t) -> 0.) where N <: AbstractArray{<:Integer}
#   func_wrapper = FunctionWrapper{Float64, Tuple{<:AbstractArray{Float64, 1}, Float64}}(func)
#   return EssentialBC{typeof(nodes), typeof(func_wrapper)}(nodes, dof, func_wrapper)
# end

# function EssentialBC(
#   nset::N,
#   dof::Int,
#   func::Function = (x, t) -> 0.
# ) where N <: NodeSet
#   return EssentialBC(nset.nodes, dof, func)
# end

# # function EssentialBC(nodes::N, dof::Int, func_string::String) where N <: AbstractArray{<:Integer}
# #   func_temp = eval(Meta.parse(func_string))
# #   func = (x, t) -> Base.invokelatest(func_temp, x, t)
# #   return EssentialBC(nodes, dof, func)
# # end

# function EssentialBC(nset::N, dof::Int, func_string::String) where N <: NodeSet
#   func_temp = eval(Meta.parse(func_string))
#   func = (x, t) -> Base.invokelatest(func_temp, x, t)
#   # return EssentialBC(nset.nodes, dof, func)
#   return EssentialBC(nset, dof, func)
# end

# # function setup_displacement_bcs(
# #   common::CthoniosCommon,
# #   exo::E,
# #   input_settings
# # ) where E <: ExodusDatabase

# #   @timeit timer(common) "Boundary Conditions" begin
# #     new_section("Boundary Conditions")
# #     @assert "boundary conditions" in keys(input_settings)
# #     @warn "Currently only supporting displacement bcs"
# #     @warn "Currently only supporting nodeset ids. Need to support names as well"

# #     bcs_settings = input_settings["boundary conditions"]
# #     @assert "displacement" in keys(bcs_settings)

# #     disp_bcs_settings = bcs_settings["displacement"]
# #     disp_bcs = EssentialBC[]
# #     for bc in disp_bcs_settings
# #       @assert "nodeset id" in keys(bc)
# #       @assert "dof" in keys(bc)
      
# #       nset = read_set(exo, NodeSet, bc["nodeset id"])
# #       dof  = bc["dof"]

# #       @info "Displacement BC"
# #       @info "  Nodeset ID      = $(nset.id)"
# #       @info "  Number of nodes = $(length(nset.nodes))"
# #       @info "  Active Dof      = $dof"

# #       # TODO add other case for tabular functions
# #       if "function" in keys(bc)
# #         func_string = bc["function"]

# #         @info "  Function        = $func_string"

# #         push!(disp_bcs, EssentialBC(nset, dof, func_string))
# #       else
# #         push!(disp_bcs, EssentialBC(nset, dof))
# #       end
# #     end
# #   end

# #   return disp_bcs
# # end
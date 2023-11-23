struct EssentialBC{V <: AbstractArray{<:Integer}, F}
  nodes::V
  dof::Int
  func::F
end

function EssentialBC(
  nset::N,
  dof::Int,
  func::Function = (x, t) -> 0.
) where N <: NodeSet
  return EssentialBC(nset.nodes, dof, func)
end

function EssentialBC(nset::N, dof::Int, func_string::String) where N <: NodeSet
  func_temp = eval(Meta.parse(func_string))
  func = (x, t) -> Base.invokelatest(func_temp, x, t)
  return EssentialBC(nset.nodes, dof, func)
end

function setup_displacement_bcs(
  common::CthoniosCommon,
  exo::E,
  input_settings
) where E <: ExodusDatabase

  @timeit timer(common) "Boundary Conditions" begin
    new_section("Boundary Conditions")
    @assert "boundary conditions" in keys(input_settings)
    @warn "Currently only supporting displacement bcs"
    @warn "Currently only supporting nodeset ids. Need to support names as well"

    bcs_settings = input_settings["boundary conditions"]
    @assert "displacement" in keys(bcs_settings)

    disp_bcs_settings = bcs_settings["displacement"]
    disp_bcs = EssentialBC[]
    for bc in disp_bcs_settings
      @assert "nodeset id" in keys(bc)
      @assert "dof" in keys(bc)
      
      nset = read_set(exo, NodeSet, bc["nodeset id"])
      dof  = bc["dof"]

      @info "Displacement BC"
      @info "  Nodeset ID      = $(nset.id)"
      @info "  Number of nodes = $(length(nset.nodes))"
      @info "  Active Dof      = $dof"

      # TODO add other case for tabular functions
      if "function" in keys(bc)
        func_string = bc["function"]

        @info "  Function        = $func_string"

        push!(disp_bcs, EssentialBC(nset, dof, func_string))
      else
        push!(disp_bcs, EssentialBC(nset, dof))
      end
    end
  end

  return disp_bcs
end
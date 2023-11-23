# make an abstract type and set up different domains for quasi-static, vs. dynamic etc.
abstract type CthoniosDomain end

struct Domain{
  Coords <: AbstractArray{<:AbstractFloat},
  Dof <: DofManager,
  Secs <: AbstractArray{<:TotalLagrangeSection},
  DispBCs <: AbstractArray{<:EssentialBC}
}
  coordinates::Coords
  dof::Dof
  sections::Secs
  disp_bcs::DispBCs
end

n_dimensions(d::Domain) = size(d.coordinates, 1)
# n_nodes_per_element(d::Domain, block_index) = size(d.connectivities[block_index], 1)
coordinates(d::Domain) = d.coordinates
# connectivity(d::Domain, block_index::Int) = d.connectivities[block_index]
# element_level_coordinates(d::Domain, block_index::Int) = @views d.coordinates[:, d.connectivities[block_index]]

# this is currently not type stable since we don't no n_dim/etc. at call time
# function element_level_coordinates_static(d::Domain, block_index::Int)
#   el_coords = element_level_coordinates(d, block_index)
#   n_dim = n_dimensions(d)
#   n_nodes = n_nodes_per_element(d, block_index)
#   reinterpret(SMatrix{n_dim, n_nodes, eltype(el_coords), n_dim * n_nodes}, vec(el_coords))
# end

# dof helper methods
create_fields(d::Domain) = create_fields(d.dof)
create_unknowns(d::Domain) = create_unknowns(d.dof)

# section helper methods
update_section_values!(d::Domain, index::Int, U::M) where M <: AbstractMatrix = 
update_section_values!(d.sections[index], U)

function update_section_values!(d::Domain, U::M) where M <: AbstractMatrix
  for n in axes(d.sections, 1)
    update_section_values!(d.sections[n], U)
  end
end

"""
"""
function Domain(common::CthoniosCommon, input_file::String)
  with_logger(common) do
    @timeit timer(common) "Domain setup" begin
      # read input file
      @timeit timer(common) "Read input file" begin
        input_settings = YAML.load_file(input_file)
      end
      
      # read from mesh

      @timeit timer(common) "Mesh" begin
        new_section("Mesh")

        @assert "mesh" in keys(input_settings)
        @assert "file name" in keys(input_settings["mesh"])
        exo = ExodusDatabase(input_settings["mesh"]["file name"], "r")
        @info "Opened ExodusDatabase"

        # read coordinates
        @timeit timer(common) "Read coordinates" begin
          @info "Reading coordinates..."
          coordinates = read_coordinates(exo)
          @info "Finished reading coordinates"
          @info "Number of dimensions = $(size(coordinates, 1))"
          @info "Number of nodes      = $(size(coordinates, 2))"
          @info ""
        end

        @timeit timer(common) "Read blocks" begin
          @warn "Currently reading all blocks by default!"
          block_ids = read_ids(exo, Block)
          blocks = Dict{Int, Block}()
          @info "Reading blocks...\n"
          for block_id in block_ids
            @info "Reading block with id $block_id"
            block = read_set(exo, Block, block_id)
            blocks[block_id] = block
            @info block
          end
          @info "Finished reading blocks"          
        end
      end

      # setup displacement bcs
      disp_bcs = setup_displacement_bcs(common, exo, input_settings)

      # TODO add stuff for traction bcs

      # setup up dof manager
      @timeit timer(common) "DofManager" begin
        new_section("DofManager")
        dof = DofManager(coordinates)
        U   = create_fields(dof)
        # update bcs for time 0.0
        @timeit timer(common) "update BCs" begin
          update_displacement_bcs!(dof, U, disp_bcs, coordinates, 0.0)
        end
      end

      # read material models
      @timeit timer(common) "Materials" begin
        new_section("Materials")
        @assert "materials" in keys(input_settings)
        material_inputs = input_settings["materials"]
        material_models = Dict{String, ConstitutiveModels.MechanicalModel}()
        material_props  = Dict{String, Vector{Float64}}()
        
        for material_name in keys(material_inputs)
          @assert "model" in keys(material_inputs[material_name])
          @assert "properties" in keys(material_inputs[material_name])
          model_name = material_inputs[material_name]["model"]
          props_dict = material_inputs[material_name]["properties"]
          @info "$material_name"
          @info "$(rpad("  model", 24)) = $model_name"
          for (key, val) in props_dict
            @info "$(rpad("  $key", 24)) = $val"
          end
          @info ""
          model_sym = Symbol(model_name)
          model_temp, props_temp = @eval $model_sym($props_dict)
          material_models[material_name] = model_temp
          material_props[material_name]  = props_temp
        end
      end

      # now for function spaces
      @timeit timer(common) "Sections" begin
        new_section("Sections")
        @assert "sections" in keys(input_settings)
        @info "Setting up sections...\n"
        section_settings = input_settings["sections"]
        sections = TotalLagrangeSection[]
        for (n, section) in enumerate(section_settings)
          @assert "block id" in keys(section)
          @assert "material" in keys(section)
          
          if "formulation" in keys(section)
            if section["formulation"] == "plane strain"
              formulation = PlaneStrain()
            else
              @assert false "Undefined section element type formulation"
            end
          else
            formulation = DefaultFormulation()
          end

          block_id = section["block id"]
          material_name = section["material"]

          @info "Section $n"
          @info "  block id = $block_id"
          @info "  material = $material_name"
          @info ""
          
          section_temp = TotalLagrangeSection(
            coordinates, 
            blocks[block_id], 
            formulation,
            material_models[material_name],
            material_props[material_name],
            common
          )
          push!(sections, section_temp)
        end
        @info "\nFinished setting up sections"
      end
    end

    # don't forget to close the ExodusDatabase
    close(exo)

    return Domain(coordinates, dof, sections, disp_bcs)
  end
end
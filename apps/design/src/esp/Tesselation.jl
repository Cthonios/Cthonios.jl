struct Tesselation{T, V}
    geometry::Geometry{T, V}
    mode::String

    function Tesselation(csm_file::String; mode::String = "python")
        return Tesselation(Geometry(csm_file), mode)
    end
    
    function Tesselation(g::Geometry{T, V}, mode::String) where {T, V}
        @assert mode == "julia" || mode == "python"
        new{T, V}(g, mode)
    end
end

function Base.close(t::Tesselation)
    close(t.geometry)
end

function tesselate!(t::Tesselation)
    base_file = split(t.geometry.csm_file, ".csm")[1]
    temp_file = base_file * "_temp.csm"
    # temp_temp_file = base_file * "_temp_temp.csm"
    tess_file = base_file * ".eto"

    if isfile(tess_file)
        # @info "Removing old tesselation file $tess_file"
        rm(tess_file)
    end

    # copy temp file
    open(temp_file, "w") do f_out
        write(f_out, "attribute capsAIM \$aflr2AIM;platoAIM\n")
        write(f_out, "attribute capsMeshLength 1.0 \n")
    
        open(t.geometry.csm_file, "r") do f_in
            for line in eachline(f_in)
                write(f_out, line, "\n")
            end
        end
    end
    # cp(temp_temp_file, temp_file, force = true)
    # rm(temp_temp_file, force = true)

    if t.mode == "python"
        _esp_python(temp_file, t.geometry.params.design_params.values)
    elseif t.mode == "julia"
        _esp_julia(temp_file)
    else
        error("Unsupported esp_mode = $(t.mode).")
    end

    _read_tesselation(t.geometry, temp_file)
end

function _esp_python(temp_file, param_vals)
    if isdir("ESP_Mesh")
        # @info "Removing old workdirs for AIMs"
        rm("ESP_Mesh", force=true, recursive=true)
    end

    script_path = joinpath(dirname(dirname(@__FILE__)), "esp", "scripts", "run_esp.bash")
    params_str = ""
    for param in param_vals
        params_str = params_str * "$param "
    end
    
    run(`$script_path $temp_file $params_str`)
    # run(`$script_path $temp_file $params_str`; stdout="/dev/null")
    # run(pipeline(`$script_path $temp_file $params_str`, stderr=devnull, stdout=devnull))
end

function _read_tesselation(g::Geometry, temp_file; verbose=false)
    counter = 0
    model = g.model
    tess_base_name = split(temp_file, ".csm")[1]
    tess_file_name = split(temp_file, ".csm")[1] * ".eto"
    exo_file_name = split(temp_file, ".csm")[1] * ".exo"

    # need to copy files over and clean things update
    out_dir = joinpath("ESP_Mesh", "Scratch")
    og_tess = joinpath(out_dir, "aflr2", "aflr2.eto")
    og_exo = joinpath(out_dir, "aflr2", "aflr2.exo")
    
    if isfile(exo_file_name)
        rm(exo_file_name, force = true)
    end
    
    if isfile(tess_file_name)
        rm(tess_file_name, force = true)
    end
    
    cp(og_exo, exo_file_name, force = true)
    cp(og_tess, tess_file_name, force = true)

    num_active_bodies = 0
    for i = 1:model.nbody
        body = unsafe_load(model.body, i + 1)
        if body.onstack == 1
            num_active_bodies = num_active_bodies + 1
        end
    end

    for i = 1:model.nbody
        body = unsafe_load(model.body, i + 1)
        if body.onstack != 1 || body.botype == OCSM_NULL_BODY
            continue
        end

        t_body = body.ebody
        if body.etess != C_NULL
            EG_deleteObject(body.etess)
        end

        # trying stuff
        curr_tess_name = tess_file_name
        if num_active_bodies > 1
            curr_tess_name = tess_base_name * "_$(counter).eto"
        end
        status = EG_loadTess(t_body, curr_tess_name, Ref(body.etess))
        _check_egads_status(status, "Failed to load tesselation file $curr_tess_name")
    end
end

struct TesselationOptions
    mode::String
end

struct ESPModel
    m_build_to::Ref{Cint}
    model::modl_T
    model_ptr::Ptr{modl_T}
end

mutable struct Geometry{
    T <: Number,
    V <: AbstractVector{T},
    D <: AbstractVector{V}
}
    csm_file::String
    dxdps::D
    params::Parameters{T, V}
    temp_csm_file::String
    tess_options::TesselationOptions
end

function Geometry(csm_file::String; mode = "python", verbose::Bool = false)
    if !verbose
        ocsmSetOutLevel(0)
    end

    # get params
    parameters = Parameters(csm_file)
    dxdps = [Vector{Float64}(undef, 0) for _ in 1:length(parameters.design_params)]
    geometry = Geometry(csm_file, dxdps, parameters, csm_file, TesselationOptions(mode))
    return geometry
end

# implementation subjcet to change
function _body_tesselation!(tess_ptr)

    points = Ref{Ptr{Cdouble}}()
    tric = Ref{Ptr{Cint}}()
    tris = Ref{Ptr{Cint}}()
    uv = Ref{Ptr{Cdouble}}()

    if tess_ptr == C_NULL
        @error "EGADS NULL OBJECT"
    end

    ref = ego()
    tess = unsafe_load(tess_ptr)

    if tess.magicnumber != MAGIC
        @error "EGADS NOT OBJ"
    end

    if tess.oclass != TESSELLATION
        @error "EGADS NOT TESS"
    end

    if tess.blind == C_NULL
        @error "EGADS NOT FOUND"
    end

    # check tesselation status
    sta = Ref{Cint}()
    n_global = Ref{Cint}()
    status = EG_statusTessBody(tess_ptr, Ref(ref), sta, n_global)
    if status != 0 || n_global[] == 0
        error("EGADS ERROR CODE: $status NGLOBAL: $(n_global[])")
    end

    btess = unsafe_load(Base.unsafe_convert(Ptr{egTessel}, tess.blind))
    n_edge, n_face = btess.nEdge, btess.nFace
    n_tri = 0

    pindex = Ref{Ptr{Cint}}()
    plen = Ref{Cint}()
    ptype = Ref{Ptr{Cint}}()
    tlen = Ref{Cint}()

    for i in 1:n_face
        status = EG_getTessFace(tess_ptr, i, plen, points, uv, ptype, pindex, tlen, tris, tric)
        if status != EGADS_SUCCESS
            @info "Face $i: EG_getTessFace status = $status (bodyTess)!\n"
        else
            n_tri = n_tri + tlen[]
        end
    end

    # get memory associated with points
    n_global = n_global[]
    xyzs = Vector{Cdouble}(undef, 3 * n_global)
    tags = Vector{verTags}(undef, n_global)

    pindex = Ref{Cint}()
    ptype = Ref{Cint}()

    # get global data
    for i = 1:n_global
        status = EG_getGlobal(tess_ptr, i, ptype, pindex, pointer(xyzs, 3 * (i - 1) + 1))
        _check_egads_status(status, "  Error: EG_getGlobal")
        tags[i] = verTags(ptype[], pindex[])
    end

    # fill up the whole triangle list -- a Face at a time
    tri = Vector{Cint}(undef, 4 * n_tri)
    

    # TODO
    n_tri = 1
    for j in 1:n_face
        pindex = Ref{Ptr{Cint}}()
        plen = Ref{Cint}()
        ptype = Ref{Ptr{Cint}}()
        tlen = Ref{Cint}()
        status = EG_getTessFace(tess_ptr, j, plen, points, uv, ptype, pindex, tlen, tris, tric)
        if status != EGADS_SUCCESS
            continue
        end

        for i in 1:tlen[]
            status = EG_localToGlobal(tess_ptr, j, unsafe_load(tris[], 3 * (i - 1) + 1), pointer(tri, 4 * (n_tri - 1) + 1))
            # _check_egads_status(status, " Face $j  $i/1 Error: EG_localToGlobal = $status (bodyTess)!\n")
            _check_egads_status(status, "  Face Error: EG_localToGlobal: coordinate y")
            status = EG_localToGlobal(tess_ptr, j, unsafe_load(tris[], 3 * (i - 1) + 2), pointer(tri, 4 * (n_tri - 1) + 2))
            # _check_egads_status(status, " Face $j  $i/2 Error: EG_localToGlobal = $status (bodyTess)!\n")
            _check_egads_status(status, "  Face Error: EG_localToGlobal: coordinate y")
            status = EG_localToGlobal(tess_ptr, j, unsafe_load(tris[], 3 * (i - 1) + 3), pointer(tri, 4 * (n_tri - 1) + 3))
            # _check_egads_status(status, " Face $j  $i/3 Error: EG_localToGlobal = $status (bodyTess)!\n")
            _check_egads_status(status, "  Face Error: EG_localToGlobal: coordinate y")
            tri[4 * (n_tri - 1) + 3 + 1] = Cint(j)
            n_tri = n_tri + 1
        end
    end

    # TODO might need to return some more vars

    n_vert = n_global
    return n_vert, n_edge, n_face, tags
end

function build_model(csm_file::String)
    # setup model
    context = Ref{ego}() # Maybe ok?
    model = Ref{Ptr{Cvoid}}()

    # open
    status = EG_open(context)
    _check_egads_status(status)
    ocsmLoad(csm_file, model)
    _check_ocsm_status(status)
    @info "Model $(csm_file) loaded"

    # setup model object
    model_ptr = Ptr{modl_T}(model[])
    model = unsafe_load(model_ptr)
    model.context = context[]

    # check file
    status = ocsmCheck(model_ptr)
    _check_ocsm_status(status)
    @info "Checks OK"

    # build geometry
    t_build_to = Cint(0)
    m_build_to = Ref{Cint}()
    n_body = Ref{Cint}()
    status = ocsmBuild(model_ptr, t_build_to, m_build_to, n_body, C_NULL)
    _check_ocsm_status(status)
    EG_deleteObject(context[])
    model.context = C_NULL
    @info "Model geometry built"

    # do some checks
    model = unsafe_load(model_ptr)
    n_body = 0
    for i in 1:model.nbody
        body = unsafe_load(model.body, i + 1)
        if body.onstack != 1 || body.botype == OCSM_NULL_BODY
            continue
        end
        n_body = n_body + 1
    end

    if n_body <= 0
        error("No bodies found.")
    end

    return ESPModel(m_build_to, model, model_ptr)
end

#TODO
function close_model(model::ESPModel)
    model_ptr = Ptr{Cvoid}(pointer_from_objref(model.model))
    ocsmFree(Ref(model_ptr))
    ocsmFree(C_NULL)
    EG_close(model.model.context)
end

# this method assumes that you have an exodus mesh already
# should probably clean this up to tesselate or check that one
# exists and is consistent
function dump_sensitivities(g::Geometry, exo_file_name)
    n_nodes = length(g.dxdps[1]) ÷ 3
    mesh_file = splitext(g.csm_file)[1] * ".exo"
    names = g.params.design_params.names
    names_all = map(x -> "d$(x)_x", names)
    append!(names_all, map(x -> "d$(x)_y", names))
    append!(names_all, map(x -> "d$(x)_z", names))

    copy_mesh(mesh_file, exo_file_name)
    ExodusDatabase(exo_file_name, "rw") do exo
        write_time(exo, 1, 0.0)
        write_names(exo, NodalVariable, names_all)

        for (n, name) in enumerate(names)
            temp_dxdp = reshape(g.dxdps[n], 3, n_nodes)
            write_values(exo, NodalVariable, 1, "d$(name)_x", temp_dxdp[1, :])
            write_values(exo, NodalVariable, 1, "d$(name)_y", temp_dxdp[2, :])
            write_values(exo, NodalVariable, 1, "d$(name)_z", temp_dxdp[3, :])
        end
    end
end

function _esp_julia(temp_file, param_vals)
    if isdir("ESP_Mesh")
        rm("ESP_Mesh", force = true, recursive = true)
    end

    prob = Problem(temp_file)
    mesh_length = get_value(get_geometry_output(prob, "MeshLength"))

    aflr2 = Analysis(prob, "aflr2AIM", "aflr2")
    set_input_value!(aflr2, "Mesh_Quiet_Flag", true)
    set_input_value!(aflr2, "Mesh_Sizing", nothing)
    set_input_value!(aflr2, "Tess_Params", [mesh_length, 1.0, 20.0])

    area_mesh = get_output_value(aflr2, "Area_Mesh")

    plato = Analysis(prob, "platoAIM", "plato")
    set_input_value!(plato, "Mesh_Morph", false)
    mesh = get_input_value(plato, "Mesh")
    link_value!(area_mesh, mesh)

    pre_analysis!(plato)
    post_analysis!(plato)

    # TODO do proper cleanup
    # close(plato)
    # close(aflr2)
    # close(prob)
    # close(plato)
end

function _esp_python(temp_file, param_vals)
    if isdir("ESP_Mesh")
        rm("ESP_Mesh", force = true, recursive = true)
    end

    script_path = joinpath(dirname(dirname(@__FILE__)), "esp", "scripts", "run_esp.bash")
    params_str = ""
    for param in param_vals
        params_str = params_str * "$param "
    end
    
    run(`$script_path $temp_file $params_str`)
end

function get_design_parameters(g::Geometry)
    return g.params.design_params.values
end

function load_tesselation(g::Geometry)
    counter = 0
    model = build_model(g.temp_csm_file)
    tess_file_name = splitext(g.csm_file)[1] * ".eto"
    exo_file_name = splitext(g.csm_file)[1] * ".exo"

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
    for i = 1:model.model.nbody
        body = unsafe_load(model.model.body, i + 1)
        if body.onstack == 1
            num_active_bodies = num_active_bodies + 1
        end
    end

    for i = 1:model.model.nbody
        body = unsafe_load(model.model.body, i + 1)
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
    return model
end

function _sensitivity(model::ESPModel, param_name; verbose::Bool = false)
    all_dxdp = Vector{Cdouble}(undef, 0)
    param_found = false
    param_index = -1

    temp_num_rows = Cint(1)
    temp_num_cols = Cint(1)
    temp_type = Cint(0)
    temp_param_name = Vector{UInt8}(undef, 129)

    for j in 1:model.model.npmtr
        temp_index = Cint(j)
        status = ocsmGetPmtr(
            model.model_ptr, temp_index, 
            Ref(temp_type), Ref(temp_num_rows), Ref(temp_num_cols), 
            temp_param_name
        )
        if status != 0
            @error "Parameter $param_name no properly set"
            error("OCSM ERROR CODE: $status")
        end
        temp_param_name_str = unsafe_string(pointer(temp_param_name))
        if param_name == temp_param_name_str
            param_found = true
            param_index = j
            break
        end
    end

    if ~param_found
        @error "Parameter $param_name not properly set"
    end

    # now actually compute the sensitivity
    ocsmSetVelD(model.model_ptr, 0, 0, 0,0.0)
    ocsmSetVelD(model.model_ptr, param_index, 1, 1, 1.0)

    # @info "Before build"
    status = ocsmBuild(model.model_ptr, Cint(0), model.m_build_to, Ref(Cint(0)), C_NULL)
    _check_ocsm_status(status)

    p_sens_ref = Ref{Ptr{Cdouble}}()

    for ibody in 1:model.model.nbody
        body = unsafe_load(model.model.body, ibody + 1)
        if body.onstack != 1
            continue
        end
        if body.botype == OCSM_NULL_BODY
            continue
        end

        n_vert, n_edge, n_face, tags = _body_tesselation!(body.etess)
        # TODO can it ever be 2D?
        dxdp = zeros(Cdouble, 3 * n_vert)

        # Nodes
        for j in 1:n_vert
            if tags[j].ptype != 0
                continue
            end
            # p_sens_ref = Ref{Ptr{Cdouble}}()
            status = ocsmGetTessVel(model.model_ptr, ibody, OCSM_NODE, tags[j].pindex, p_sens_ref)
            # _check_egads_status(status, "VERTEX: $j")
            _check_egads_status(status, "VERTEX: ")
            p_sens = p_sens_ref[]
            dxdp[3 * (j - 1) + 1] = unsafe_load(p_sens, 1)
            dxdp[3 * (j - 1) + 2] = unsafe_load(p_sens, 2)
            dxdp[3 * (j - 1) + 3] = unsafe_load(p_sens, 3)
        end

        # edges
        for j in 1:n_edge
            # p_sens_ref = Ref{Ptr{Cdouble}}()
            status = ocsmGetTessVel(model.model_ptr, ibody, OCSM_EDGE, j, p_sens_ref)
            # _check_egads_status(status, "EDGE: $j")
            _check_egads_status(status, "EDGE: ")
            p_sens = p_sens_ref[]

            for k in 1:n_vert
                if tags[k].ptype > 0 && tags[k].pindex == j
                    temp_index = tags[k].ptype - 1
                    dxdp[3 * (k - 1) + 1] = unsafe_load(p_sens, 3 * temp_index + 1)
                    dxdp[3 * (k - 1) + 2] = unsafe_load(p_sens, 3 * temp_index + 2)
                    dxdp[3 * (k - 1) + 3] = unsafe_load(p_sens, 3 * temp_index + 3)
                end
            end
        end
        
        # faces
        for j in 1:n_face
            # p_sens_ref = Ref{Ptr{Cdouble}}()
            status = ocsmGetTessVel(model.model_ptr, ibody, OCSM_FACE, j, p_sens_ref)
            # _check_egads_status(status, "FACE: $j")
            _check_egads_status(status, "FACE: ")
            p_sens = p_sens_ref[]

            for k in 1:n_vert
                if tags[k].ptype < 0 && tags[k].pindex == j
                    temp_index = -tags[k].ptype - 1
                    dxdp[3 * (k - 1) + 1] = unsafe_load(p_sens, 3 * temp_index + 1)
                    dxdp[3 * (k - 1) + 2] = unsafe_load(p_sens, 3 * temp_index + 2)
                    dxdp[3 * (k - 1) + 3] = unsafe_load(p_sens, 3 * temp_index + 3)
                end
            end
        end
        append!(all_dxdp, dxdp)
    end
    return all_dxdp
end

function tesselate(g::Geometry; verbose::Bool = false)
    mktempdir() do path
        temp_file = joinpath(path, basename(g.temp_csm_file))
        open(temp_file, "w") do f_out
            write(f_out, "attribute capsAIM \$aflr2AIM;platoAIM\n")
            write(f_out, "attribute capsMeshLength 1.0 \n")
    
            open(g.temp_csm_file, "r") do f_in
                for line in eachline(f_in)
                    write(f_out, line, "\n")
                end
            end
        end

        if verbose
            open(temp_file, "r") do f
                lines = readlines(f)
                for line in lines
                    println(line)
                end
            end
        end

        if g.tess_options.mode == "python"
            _esp_python(temp_file, get_design_parameters(g))
        elseif g.tess_options.mode == "julia"
            _esp_julia(temp_file, get_design_parameters(g))
        else
            error("Unsupported esp_mode = $(g.tess_options.mode).")
        end
    end
end

# TODO make an option to save in a sepcific dir name maybe?
function update!(g::Geometry, p; verbose::Bool = false)
    @assert length(p) == length(g.params.design_params)
    g.params.design_params.values .= p
    mktempdir() do path
        counter = 1
        temp_file = joinpath(path, "temp.csm")
        open(temp_file, "w") do outf
            for line in eachline(g.csm_file)
                stripped = strip(line)
                if startswith(stripped, "despmtr")
                    parts = split(stripped)
                    pname = String(parts[2])
                    @assert pname in g.params.design_params.names
                    index = findfirst(p -> p == pname, g.params.design_params.names)
                    parts[3] = @sprintf("%.16g", g.params.design_params.values[index])
                    counter = counter + 1
                    println(outf, join(parts, " "))
                else
                    println(outf, line)
                end
            end
        end

        if verbose
            for line in eachline(temp_file)
                println(line)
            end
        end

        # set this for now, update methods later
        # to make geometry immutable
        g.temp_csm_file = temp_file

        tesselate(g)
        model = load_tesselation(g)

        _update_sensitivities!(g, model; verbose = verbose)

        # TODO do proper cleanup
        # close_model(model)
    end
end

function _update_sensitivities!(g::Geometry, model::ESPModel; verbose::Bool = false)
    for (n, param_name) in enumerate(g.params.design_params.names)
        dxdp_temp = _sensitivity(model, param_name; verbose = verbose)
        g.dxdps[n] = dxdp_temp

        if verbose
            @info "CoordinateSensitivity for parameter $param_name calculated"
        end
    end
end

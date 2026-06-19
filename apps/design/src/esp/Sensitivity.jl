struct CoordinateSensitivity{T, A <: AbstractArray{<:AbstractArray{T, 1}, 1}}
    dxdps::A
    tesselation::Tesselation
end

function CoordinateSensitivity(t::Tesselation; verbose=false)
    # TODO currently reading first param sensitivity once to get sizes
    temp = _sensitivity(
        t.geometry.m_build_to, 
        t.geometry.model, 
        t.geometry.model_ptr, 
        # t.geometry.design_parameter_names[1]
        t.geometry.params.design_params.names[1]
    )
    dxdps = Vector{Vector{Cdouble}}(undef, length(t.geometry.params.design_params.names))
    _sensitivities_init!(dxdps, t)

    n_params = length(t.geometry.params.design_params.names)    
    n_nodes = length(temp) ÷ 3

    @assert length(dxdps) == n_params
    for n in axes(dxdps, 1)
        @assert length(dxdps[n]) == 3 * n_nodes
    end
    return CoordinateSensitivity{Float64, typeof(dxdps)}(dxdps, t)
end

function CoordinateSensitivity(csm_file; verbose=false)
    tesselation = Tesselation(csm_file; mode="python")
    return CoordinateSensitivity(tesselation; verbose=verbose)
end

function Base.close(s::CoordinateSensitivity)
    close(s.tesselation.geometry)
end

function Base.show(io::IO, s::CoordinateSensitivity)
    println(io, "Coordinate sensitivty:")
    println(io, "  Number of variables    = $(length(s.dxdps))")
    println(io, "  Sensitivty vector size = $(size(s.dxdps[1]))")
    show(io, s.tesselation.geometry.params; pad = "  ")
end

function dump_sensitivities(s::CoordinateSensitivity, exo_file_name)
    
    n_params = length(s.tesselation.geometry.design_parameter_names)
    n_nodes = length(s.dxdps[1]) ÷ 3

    # dxdp = reshape(s.dxdp, n_params, 3, n_nodes)
    mesh_file = split(s.tesselation.geometry.csm_file, ".")[1] * "_temp.exo"
    names = s.tesselation.geometry.design_parameter_names
    names_all = map(x -> "d$(x)_x", names)
    append!(names_all, map(x -> "d$(x)_y", names))
    append!(names_all, map(x -> "d$(x)_z", names))

    copy_mesh(mesh_file, exo_file_name)
    exo = ExodusDatabase(exo_file_name, "rw")
    write_time(exo, 1, 0.0)
    write_names(exo, NodalVariable, names_all)

    for (n, name) in enumerate(names)
        temp_dxdp = reshape(s.dxdps[n], 3, n_nodes)
        write_values(exo, NodalVariable, 1, "d$(name)_x", temp_dxdp[1, :])
        write_values(exo, NodalVariable, 1, "d$(name)_y", temp_dxdp[2, :])
        write_values(exo, NodalVariable, 1, "d$(name)_z", temp_dxdp[3, :])
    end
    close(exo)
end

function update_design_parameters!(s::CoordinateSensitivity, param_vals)
    update_design_parameters!(s.tesselation, param_vals)
    _sensitivities_update!(s.dxdps, s.tesselation)
    return nothing
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
        # @time _check_egads_status(status, " Error: EG_getGlobal $i status = $status (bodyTessellation)!\n")
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

# Need to two methods since the memory isn't set
function _sensitivities_init!(dxdps, t::Tesselation; verbose=false)
    for (n, param_name) in enumerate(t.geometry.params.design_params.names)
        dxdp_temp = _sensitivity(
            t.geometry.m_build_to, 
            t.geometry.model, 
            t.geometry.model_ptr, 
            param_name; verbose=verbose
        )
        dxdps[n] = dxdp_temp

        if verbose
            @info "CoordinateSensitivity for parameter $param_name calculated"
        end
    end
end

function _sensitivities_update!(dxdps, t::Tesselation; verbose=false)
    for (n, param_name) in enumerate(t.geometry.design_parameter_names)
        dxdp_temp = _sensitivity(
            t.geometry.m_build_to, 
            t.geometry.model, 
            t.geometry.model_ptr, 
            param_name; verbose=verbose
        )
        resize!(dxdps[n], length(dxdp_temp))
        copyto!(dxdps[n], dxdp_temp)

        if verbose
            @info "CoordinateSensitivity for parameter $param_name calculated"
        end
    end
end

function _sensitivity(m_build_to, model, model_ptr, param_name; verbose=false)
    all_dxdp = Vector{Cdouble}(undef, 0)
    param_found = false
    param_index = -1

    temp_num_rows = Cint(1)
    temp_num_cols = Cint(1)
    temp_type = Cint(0)
    temp_param_name = Vector{UInt8}(undef, 129)

    for j in 1:model.npmtr
        temp_index = Cint(j)
        # temp_num_rows = Cint(1)
        # temp_num_cols = Cint(1)
        # temp_type = Cint(0)
        # temp_param_name = Vector{UInt8}(undef, 129)
        status = ocsmGetPmtr(
            model_ptr, temp_index, 
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
    ocsmSetVelD(model_ptr, 0, 0, 0,0.0)
    ocsmSetVelD(model_ptr, param_index, 1, 1, 1.0)

    # @info "Before build"
    status = ocsmBuild(model_ptr, Cint(0), m_build_to, Ref(Cint(0)), C_NULL)
    _check_ocsm_status(status)

    p_sens_ref = Ref{Ptr{Cdouble}}()

    for ibody in 1:model.nbody
        body = unsafe_load(model.body, ibody + 1)
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
            status = ocsmGetTessVel(model_ptr, ibody, OCSM_NODE, tags[j].pindex, p_sens_ref)
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
            status = ocsmGetTessVel(model_ptr, ibody, OCSM_EDGE, j, p_sens_ref)
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
            status = ocsmGetTessVel(model_ptr, ibody, OCSM_FACE, j, p_sens_ref)
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

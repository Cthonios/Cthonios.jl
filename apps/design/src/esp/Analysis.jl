struct CapsInputInfo
    vtype::capsvType
    nrow::Int
    ncol::Int
    units::Union{Nothing, String}
end

function _caps_info(obj::capsObj)
    name    = Ref{Ptr{Cchar}}()
    type    = Ref{capsoType}()
    subtype = Ref{capssType}()
    link    = Ref{capsObj}()
    parent  = Ref{capsObj}()
    last    = Ref{capsOwn}()
    ierr = caps_info(obj, name, type, subtype, link, parent, last)

    if ierr != 0
        @error "CAPS INFO ERROR CODE: $ierr"
    end

    return unsafe_string(name[]), type[], subtype[],
        link[], parent[], last[]
end

function _caps_problem(csm_file, pr_name)
    csm_file_c = Vector{Cchar}(collect(csm_file))
    n_err = Ref{Cint}()
    caps = Ref{capsObj}()
    errs = Ref{Ptr{capsErrs}}()
    status = caps_open(
        pr_name, C_NULL, oFileName, 
        pointer(csm_file_c), Cint(0), caps,
        n_err, errs
    )
    if status != 0
        error("CAPS ERROR CODE: $status")
    end
    return caps
end

function Base.close(problem::capsObj)
    status = caps_close(problem, Cint(0), C_NULL)
    if status != 0
        @error "Error inc lose caps obj: $status"
    end
end

function convert_to_caps(info::CapsInputInfo, value)
    if info.vtype == Boolean
        value isa Bool || error("Expected Bool, got $(typeof(value))")
        storage = Ref{Cint}(value ? 1 : 0)
        return (
            info.vtype,
            1,
            1,
            Base.unsafe_convert(Ptr{Cvoid}, storage),
            storage,
        )
    elseif info.vtype == Double
        if value isa Real
            storage = Ref{Cdouble}(Float64(value))
            return (
                info.vtype,
                1,
                1,
                Base.unsafe_convert(Ptr{Cvoid}, storage),
                storage,
            )
        elseif value isa AbstractVector{<:Real}
            storage = Vector{Cdouble}(Float64.(value))
            return (
                info.vtype,
                length(storage),
                1,
                # pointer(storage),
                # Ref(storage),
                pointer(storage),
                # storage,
                storage,
            )
        else
            error("Cannot convert $(typeof(value)) to CAPS Double")
        end
    elseif info.vtype == mString
        value isa String || error("Expected String")
        storage = Base.cconvert(Cstring, value)
        return (
            info.vtype,
            1,
            1,
            Ptr{Cvoid}(storage),
            storage,
        )
    elseif info.vtype == Tuple
        value isa Nothing || error("Expected Nothing. Need to support tuples")
        storage = Ptr{C_NULL}(C_NULL)
        return (
            info.vtype,
            1,
            1,
            storage,
            storage
        )
    else
        error("Conversion for CAPS type $(info.vtype) not implemented")
    end
end

function input_info(obj::capsObj)
    vtype = Ref{capsvType}()
    nrow  = Ref{Cint}()
    ncol  = Ref{Cint}()

    data    = Ref{Ptr{Cvoid}}()
    partial = Ref{Ptr{Cint}}()
    units   = Ref{Ptr{Cchar}}()

    n_err = Ref{Cint}()
    errs  = Ref{Ptr{capsErrs}}()

    status = caps_getValue(
        obj,
        vtype,
        nrow,
        ncol,
        data,
        partial,
        units,
        n_err,
        errs,
    )

    status == 0 || error("CAPS ERROR CODE: $status")

    return CapsInputInfo(
        vtype[],
        Int(nrow[]),
        Int(ncol[]),
        units[] == C_NULL ? nothing : unsafe_string(units[])
    )
end

function get_value(obj::capsObj)
    vtype = Ref{capsvType}()
    nrow  = Ref{Cint}()
    ncol  = Ref{Cint}()

    data    = Ref{Ptr{Cvoid}}()
    partial = Ref{Ptr{Cint}}()
    units   = Ref{Ptr{Cchar}}()

    n_err = Ref{Cint}()
    errs  = Ref{Ptr{capsErrs}}()

    status = caps_getValue(
        obj, vtype, nrow, ncol,
        data, partial, units,
        n_err, errs
    )
    status == 0 || error("CAPS ERROR CODE: $status")

    if vtype[] == Boolean
        data = unsafe_load(Base.unsafe_convert(Ptr{Bool}, data[]))
    elseif vtype[] == Integer
        if ncol[] > 1
            @assert false "Probably not supported"
        else
            data = unsafe_wrap(Vector{Int32}, Base.unsafe_convert(Ptr{Int32}, data[]), nrow[])
        end
        if length(data) == 1
            data = data[1]
        end
    # is the deriv right?
    elseif vtype[] == Double || vtype[] == DoubleDeriv
        if ncol[] > 1
            @assert false "Probably not supported"
        else
            data = unsafe_wrap(Vector{Float64}, Base.unsafe_convert(Ptr{Float64}, data[]), nrow[])
        end
        if length(data) == 1
            data = data[1]
        end
    elseif vtype[] == mString
        if data[] == Ptr{C_NULL}(C_NULL)
            data = ""
        else
            data = unsafe_string(Base.unsafe_convert(Ptr{Cchar}, data[]))
        end
    elseif vtype[] == Tuple
        data = data[]
        if data == Ptr{C_NULL}(C_NULL)
            data = nothing
        else
            @assert false "Unsupported tuple type"
        end
    else
        @assert false "Not implemented yet capsvType $(vtype[]) with data = $(data)"
    end
    return data
end

abstract type AbstractCapsWrapper end

function Base.close(an::AbstractCapsWrapper)
    close(an.obj[])
end

function Base.size(an::AbstractCapsWrapper, otype, stype)
    size = Ref{Cint}()
    n_err = Ref{Cint}()
    errs = Ref{Ptr{capsErrs}}()
    status = caps_size(an.obj[], otype, stype, size, n_err, errs)
    if status != 0
        error("CAPS ERROR CODE: $status")
    end
    size = size[]
    return size
end

function children(an::AbstractCapsWrapper, otype, stype)
    val_objs = []
    for i in 1:size(an, otype, stype)
        child = Ref{capsObj}()
        status = caps_childByIndex(an.obj[], otype, stype, Cint(i), child)
        if status != 0
            error("CAPS ERROR CODE: $status")
        end
        push!(val_objs, child[])
    end
    return val_objs
end

struct Problem <: AbstractCapsWrapper
    obj::Ref{capsObj}
end

function Problem(csm_file::String; name = "ESP_Mesh")
    obj = _caps_problem(csm_file, name)
    return Problem(obj)
end

function get_geometry_output(prob::Problem, name::String)
    obj = Ref{capsObj}()

    n_err = Ref{Cint}()
    errs  = Ref{Ptr{capsErrs}}()

    status = caps_childByName(
        prob.obj[],
        VALUE,
        GEOMETRYOUT,
        name,
        obj,
        n_err,
        errs,
    )

    status == 0 || error("CAPS ERROR CODE: $status")

    return obj[]
end

struct Analysis <: AbstractCapsWrapper
    problem::Problem
    obj::Ref{capsObj}
end

function Analysis(prob::Problem, aim_name, name; verbose=false)
    analysis_obj = Ref{capsObj}()

    intent   = C_NULL
    unit_sys = C_NULL
    execute  = Ref{Cint}(1)

    n_err = Ref{Cint}()
    errs  = Ref{Ptr{capsErrs}}()

    status = caps_makeAnalysis(
        prob.obj[],
        aim_name,
        name,
        unit_sys,
        intent,
        execute,
        analysis_obj,
        n_err,
        errs,
    )

    status == 0 || error("CAPS ERROR CODE: $status")

    return Analysis(prob, analysis_obj)
end

function dirty_analyses(an)
    na_obj = Ref{Cint}()
    a_objs = Ref{Ptr{capsObj}}()
    status = caps_dirtyAnalysis(an.obj[], na_obj, a_objs)
    if status != 0
        error("CAPS ERROR CODE: $status")
    end
    
    @info "Number of dirty analyses = $(na_obj[])"
    for a in unsafe_wrap(Vector{capsObj}, a_objs[], na_obj[])
        @show unsafe_string(unsafe_load(a).name)
    end
end

function get_input_value(an::Analysis, name::String)
    ch = children(an, VALUE, ANALYSISIN)
    ch_objs = map(Base.unsafe_load, ch)
    names = map(x -> unsafe_string(x.name), ch_objs)
    indices = findall(==(name), names)
    @assert length(indices) == 1
    index = indices[1]
    return ch[index]
end

function get_output_value(an::Analysis, name::String)
    ch = children(an, VALUE, ANALYSISOUT)
    ch_objs = map(Base.unsafe_load, ch)
    names = map(x -> unsafe_string(x.name), ch_objs)
    indices = findall(==(name), names)
    @assert length(indices) == 1
    index = indices[1]
    return ch[index]
end

function link_value!(source::capsObj, destination::capsObj)
    n_err = Ref{Cint}()
    errs = Ref{Ptr{capsErrs}}()
    tmethod = Copy
    status = caps_linkValue(
        source, tmethod, destination,
        n_err, errs
    )
    if status != 0
        error("CAPS ERROR CODE: $status")
    end
end

function post_analysis!(an::Analysis)
    n_err = Ref{Cint}()
    errs = Ref{Ptr{capsErrs}}()
    status = caps_postAnalysis(an.obj[], n_err, errs)
    if status != 0
        error("CAPS ERROR CODE: $status")
    end
    return nothing
end

function pre_analysis!(an::Analysis)
    n_err = Ref{Cint}()
    errs = Ref{Ptr{capsErrs}}()
    status = caps_preAnalysis(an.obj[], n_err, errs)
    if status != 0
        error("CAPS ERROR CODE: $status")
    end
    return nothing
end

function set_input_value!(an::Analysis, name::String, value)
    obj = get_input_value(an, name)
    info = input_info(obj)
    @info "Setting CAPS input" name info.vtype info.nrow info.ncol
    vtype, nrow, ncol, data_ptr, storage = convert_to_caps(info, value)
    # TODO
    units = C_NULL
    n_err = Ref{Cint}()
    errs  = Ref{Ptr{capsErrs}}()

    GC.@preserve storage begin
        status = caps_setValue(
            obj, vtype, Cint(nrow), Cint(ncol),
            data_ptr, C_NULL, # never use partial updates
            units,
            n_err, errs,
        )
        status == 0 || error("CAPS ERROR CODE: $status")
    end

    return nothing
end

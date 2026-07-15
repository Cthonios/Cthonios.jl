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

# function _caps_problem(csm_file, pr_name)
#     # csm_file_c = Vector{Cchar}(collect(csm_file))
#     csm_file_c = push!(Vector{Cchar}(collect(csm_file)), 0)
#     @show csm_file_c
#     n_err = Ref{Cint}()
#     caps = Ref{capsObj}()
#     # caps = capsObj()
#     errs = Ref{Ptr{capsErrs}}()
#     # # GC.@preserve csm_file_c begin
#     #     status = caps_open(
#     #         pr_name, C_NULL, oFileName, 
#     #         csm_file_c, Cint(0), caps,
#     #         n_err, errs
#     #     )
#     # # end
#     csm_file_conv = Base.cconvert(Cstring, csm_file)  # GC root; guarantees NUL-terminated buffer
#     # GC.@preserve csm_file_conv begin
#         filename_ptr = Ptr{Cvoid}(Base.unsafe_convert(Cstring, csm_file_conv))
#         status = caps_open(
#             pr_name, C_NULL, oFileName,
#             filename_ptr, Cint(0), caps,
#             n_err, errs
#         )
#     # end
#     @show caps[]
#     if status != 0
#         error("CAPS ERROR CODE: $status")
#     end
#     return caps[]
# end

# function Base.close(problem::capsObj)
#     status = caps_close(problem, Cint(1), C_NULL)
#     if status != 0
#         @error "Error in close caps obj: $status"
#     end
# end
# const _CAPS_STRING_KEEPALIVE = Any[]  # module-level, never collected

# function _caps_problem(csm_file::String, pr_name::String)
#     n_err = Ref{Cint}()
#     caps  = Ref{capsObj}()
#     errs  = Ref{Ptr{capsErrs}}()

#     csm_file_conv = Base.cconvert(Cstring, csm_file)  # NUL-terminated, GC root
#     push!(_CAPS_STRING_KEEPALIVE, csm_file_conv)   # never let GC touch this
#     # GC.@preserve csm_file_conv begin
#         filename_ptr = Ptr{Cvoid}(Base.unsafe_convert(Cstring, csm_file_conv))
#         status = caps_open(
#             pr_name, C_NULL, oFileName,
#             filename_ptr, Cint(0), caps,
#             n_err, errs
#         )
#     # end
#     status == 0 || error("CAPS ERROR CODE: $status")

#     return caps[]   # unwrap: capsObj, not Ref{capsObj}
# end

# function _caps_problem(csm_file::String, pr_name::String)
#     n_err = Ref{Cint}()
#     caps  = Ref{capsObj}()
#     errs  = Ref{Ptr{capsErrs}}()

#     csm_conv = Base.cconvert(Cstring, csm_file)
#     name_conv = Base.cconvert(Cstring, pr_name)

#     filename_ptr = Ptr{Cvoid}(Base.unsafe_convert(Cstring, csm_conv))
#     # name_ptr     = Base.unsafe_convert(Cstring, name_conv)
#     name_ptr     = pointer(name_conv)
#     status = GC.@preserve csm_conv name_conv begin
#         caps_open(name_ptr, C_NULL, oFileName, filename_ptr, Cint(1), caps, n_err, errs)
#     end

#     status == 0 || error("CAPS ERROR CODE: $status")

#     return caps[], Any[csm_conv, name_conv]   # hand the roots back up
# end
# function _caps_problem(csm_file::String, pr_name::String)
#     # Save current signal dispositions before calling into CAPS
#     old_sigint  = ccall(:signal, Ptr{Cvoid}, (Cint, Ptr{Cvoid}), 2,  C_NULL)  # SIGINT = 2 on Linux
#     old_sigsegv = ccall(:signal, Ptr{Cvoid}, (Cint, Ptr{Cvoid}), 11, C_NULL) # SIGSEGV = 11
#     old_sigabrt = ccall(:signal, Ptr{Cvoid}, (Cint, Ptr{Cvoid}), 6,  C_NULL) # SIGABRT = 6

#     n_err = Ref{Cint}()
#     caps  = Ref{capsObj}()
#     errs  = Ref{Ptr{capsErrs}}()

#     csm_conv  = Base.cconvert(Cstring, csm_file)
#     # name_conv = Base.cconvert(Cstring, pr_name)
#     filename_ptr = Ptr{Cvoid}(Base.unsafe_convert(Cstring, csm_conv))
#     # name_ptr     = Base.unsafe_convert(Cstring, name_conv)

#     status = GC.@preserve csm_conv begin
#         caps_open(pr_name, C_NULL, oFileName, filename_ptr, Cint(1), caps, n_err, errs)
#     end

#     # Restore Julia's original handlers immediately — don't let CAPS's handlers persist
#     ccall(:signal, Ptr{Cvoid}, (Cint, Ptr{Cvoid}), 2,  old_sigint)
#     ccall(:signal, Ptr{Cvoid}, (Cint, Ptr{Cvoid}), 11, old_sigsegv)
#     ccall(:signal, Ptr{Cvoid}, (Cint, Ptr{Cvoid}), 6,  old_sigabrt)

#     status == 0 || error("CAPS ERROR CODE: $status")
#     return caps[]#, Any[csm_conv, name_conv]
# end
# Minimal struct matching Linux x86_64 struct sigaction layout well enough
# to read/restore opaquely (we never interpret the fields ourselves).
const SIGACTION_SIZE = 152  # generous upper bound; real struct is smaller but padding varies

function _get_sigaction(signum::Cint)
    old = Vector{UInt8}(undef, SIGACTION_SIZE)
    ccall(:sigaction, Cint, (Cint, Ptr{Cvoid}, Ptr{UInt8}), signum, C_NULL, old)
    return old
end

function _set_sigaction(signum::Cint, act::Vector{UInt8})
    ccall(:sigaction, Cint, (Cint, Ptr{UInt8}, Ptr{Cvoid}), signum, act, C_NULL)
end

function _caps_problem(csm_file::String, pr_name::String)
    SIGSEGV = Cint(11); SIGABRT = Cint(6); SIGINT = Cint(2); SIGBUS = Cint(7); SIGHUP = Cint(1)
    saved = Dict(s => _get_sigaction(s) for s in (SIGSEGV, SIGABRT, SIGINT, SIGBUS, SIGHUP))

    n_err = Ref{Cint}(); caps = Ref{capsObj}(); errs = Ref{Ptr{capsErrs}}()
    csm_conv  = Base.cconvert(Cstring, csm_file)
    # name_conv = Base.cconvert(Cstring, pr_name)
    filename_ptr = Ptr{Cvoid}(Base.unsafe_convert(Cstring, csm_conv))
    # name_ptr     = Base.unsafe_convert(Cstring, name_conv)

    status = GC.@preserve csm_conv begin
        caps_open(pr_name, C_NULL, oFileName, filename_ptr, Cint(0), caps, n_err, errs)
    end

    # restore Julia's exact original handlers (flags, mask, and all)
    for (s, act) in saved
        _set_sigaction(s, act)
    end

    status == 0 || error("CAPS ERROR CODE: $status")
    return caps[]#, Any[csm_conv, name_conv]
end

function Base.close(problem::capsObj)
    status = caps_close(problem, Cint(0), C_NULL)
    status == 0 || @error "Error in close caps obj: $status"
    return nothing
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
        storage = value  # keep the String alive as the GC root
        data_ptr = Ptr{Cvoid}(Base.unsafe_convert(Cstring, storage))
        return (info.vtype, 1, 1, data_ptr, storage)
        # storage = Base.cconvert(Cstring, value)
        # return (
        #     info.vtype,
        #     1,
        #     1,
        #     Ptr{Cvoid}(storage),
        #     storage,
        # )
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
            new_data = Vector{Int32}(undef, nrow[])
            temp_data = unsafe_wrap(Vector{Int32}, Base.unsafe_convert(Ptr{Int32}, data[]), nrow[])
            for n in eachindex(new_data)
                new_data[n] = temp_data[n]
            end
            data = new_data
        end
        if length(data) == 1
            data = data[1]
        end
    # is the deriv right?
    elseif vtype[] == Double || vtype[] == DoubleDeriv
        if ncol[] > 1
            @assert false "Probably not supported"
        else
            new_data = Vector{Float64}(undef, nrow[])
            temp_data = unsafe_wrap(Vector{Float64}, Base.unsafe_convert(Ptr{Float64}, data[]), nrow[])
            for n in eachindex(new_data)
                new_data[n] = temp_data[n]
            end
            data = new_data
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
    close(an.obj)
end

function Base.size(an::AbstractCapsWrapper, otype, stype)
    size = Ref{Cint}()
    n_err = Ref{Cint}()
    errs = Ref{Ptr{capsErrs}}()
    status = caps_size(an.obj, otype, stype, size, n_err, errs)
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
        status = caps_childByIndex(an.obj, otype, stype, Cint(i), child)
        if status != 0
            error("CAPS ERROR CODE: $status")
        end
        push!(val_objs, child[])
    end
    return val_objs
end

struct Problem <: AbstractCapsWrapper
    obj::capsObj
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
        prob.obj,
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
    obj::capsObj
end

function Analysis(prob::Problem, aim_name, name; verbose=false)
    analysis_obj = Ref{capsObj}()

    intent   = C_NULL
    unit_sys = C_NULL
    execute  = Ref{Cint}(1)

    n_err = Ref{Cint}()
    errs  = Ref{Ptr{capsErrs}}()

    status = caps_makeAnalysis(
        prob.obj,
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

    return Analysis(prob, analysis_obj[])
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
    status = caps_postAnalysis(an.obj, n_err, errs)
    if status != 0
        error("CAPS ERROR CODE: $status")
    end
    return nothing
end

function pre_analysis!(an::Analysis)
    n_err = Ref{Cint}()
    errs = Ref{Ptr{capsErrs}}()
    status = caps_preAnalysis(an.obj, n_err, errs)
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

    # GC.@preserve storage begin
        status = caps_setValue(
            obj, vtype, Cint(nrow), Cint(ncol),
            data_ptr, C_NULL, # never use partial updates
            units,
            n_err, errs,
        )
        status == 0 || error("CAPS ERROR CODE: $status")
    # end

    return nothing
end

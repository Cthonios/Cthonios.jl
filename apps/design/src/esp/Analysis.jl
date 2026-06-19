# struct CapsProblem{Obj}
#     obj::Obj
# end

function _caps_problem(csm_file, pr_name)
    ph_name = C_NULL

    if isdir(pr_name)
        @info "Removing old AIM directory"
        rm(pr_name, force=true, recursive=true)
    end

    flag = Cint(oFileName)
    pptr = Ptr{Cvoid}(pointer(Vector{UInt8}(csm_file)))
    out_level = Cint(0)
    problem_obj = Ref{capsObj}()
    n_err = Ref{Cint}()
    errs = Ref{Ptr{capsErrs}}()
    status = caps_open(pr_name, ph_name, flag, pptr, out_level, problem_obj, n_err, errs)
    if status != 0
        error("CAPS ERROR CODE: $status")
    end

    return problem_obj
end

struct Analysis
    obj::Ref{capsObj}
end

function Analysis(problem_obj, aim_name, name; verbose=false)
    # problem_obj = prob.obj
    intent = Ptr{Int8}(C_NULL)
    unit_sys = Ptr{Int8}(C_NULL)
    execute = Cint(1)
    analysis_obj = Ref{capsObj}()
    n_err = Ref(Cint(0))
    errs = Ref{Ptr{capsErrs}}()
    status = caps_makeAnalysis(
        problem_obj[], aim_name, name, unit_sys, intent,
        Ref(execute), analysis_obj,
        n_err, errs
    )
    if status != 0
        error("CAPS ERROR CODE: $status")
    end

    # @show unsafe_load(analysis[]).info
    # get info # TODO make optional
    if verbose
        name = Ref{Ptr{Cchar}}()
        otype = Ref{capsoType}()
        stype = Ref{capssType}()
        link = Ref{capsObj}()
        parent = Ref{capsObj}()
        last = Ref{capsOwn}()
        status = caps_info(analysis_obj[], name, otype, stype, link, parent, last)
        if status != 0
            error("CAPS ERROR CODE: $status")
        end

        str = "AIM info:\n"
        str = str * "  name = $(Base.unsafe_string(name[]))"
        @info str
    end

    return Analysis(analysis_obj)
end 

function children(an::Analysis, otype, stype)
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

function _convert_data(data)
    type = typeof(data)
    vtype = get_value_type(data)
    n_row, n_col = get_value_size(data)

    # TODO
    # if vtype == Boolean || vtype == Integers || vtype == DoubleDeriv
    if type == Bool || type == Integer || type <: AbstractFloat
        partial = [Cint(1)]
        val = [data]
    elseif type <: AbstractArray
        partial = Vector{Cint}(undef, length(data))
        partial .= Cint(1)
        val = data
    else
        @warn "$vtype is still experimental"
        partial = [Cint(1)]
        val = [data]
        # @assert false "vtype unsupported"
    end
    return vtype, n_row, n_col, val, partial
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

function get_value_size(val)
    type = typeof(val)
    n_row = length(val)
    if type == Bool
        n_col = 1
    elseif type <: AbstractArray
        # data = [data]
        # n_col = length(val)
        n_col = length(val[1])
    elseif type == Ptr{Nothing}
        n_col = 1
    end

    return n_row, n_col
end

function get_value_type(val)
    type = typeof(val)
    if type <: AbstractArray
        vtype = get_value_type(val[1])
        for i in axes(val, 1)
            if vtype != get_value_type(val[i])
                error("Types not equivalent")
            end
        end
    elseif type == Bool
        vtype = Boolean
    # elseif type <: Integer
        # vtype = Integers
    elseif type <: Real
        vtype = Double
    elseif type <: String
        vtype = mString
    elseif type == Ptr{Nothing}
        # @error "DO me"
        vtype = Integers
    end
    return vtype
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

function Base.size(an::Analysis, otype, stype)
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

function set_input_value!(an::Analysis, name, data)
    val = get_input_value(an, name)

    units = Ptr{Int8}(C_NULL)
    
    # convert data
    type = typeof(data)
    if type != Bool && type != Ptr{Nothing}
        if type <: AbstractArray
            # do nothing
        else
            error("Data type $type not supported")
        end
    end

    vtype, n_row, n_col, pdata, partial = _convert_data(data)

    # errors
    n_err = Ref(Cint(0))
    errs = Ref{Ptr{capsErrs}}()

    status = caps_setValue(
        val, vtype, 
        n_row, n_col,
        pdata, partial, units, 
        n_err, errs
    )
    if status != 0
        error("CAPS ERROR CODE: $status")
    end

    return nothing
end

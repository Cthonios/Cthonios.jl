struct Geometry{T, V}
    csm_file::String
    m_build_to::Ref{Cint}
    model::modl_T
    model_ptr::Ptr{modl_T}
    params::Parameters{T, V}
end

function Geometry(csm_file::String; verbose::Bool = false)
    if !verbose
        ocsmSetOutLevel(0)
    end

    # get params
    parameters = Parameters(csm_file)

    # setup model
    context = Ref{ego}() # Maybe ok?
    model = Ref{Ptr{Cvoid}}()

    # open
    status = EG_open(context)
    _check_egads_status(status)
    ocsmLoad(csm_file, model)
    _check_ocsm_status(status)
    @info "Model $csm_file loaded"

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

    geometry = Geometry(csm_file, m_build_to, model, model_ptr, parameters)
    return geometry
end

function Base.close(g::Geometry)
    model_ptr = Ptr{Cvoid}(pointer_from_objref(g.model))
    ocsmFree(Ref(model_ptr))
    ocsmFree(C_NULL)
    EG_close(g.model.context)
end

function update!(g::Geometry, p::Parameters)
    mktemp() do path, io
        counter = 1
        temp_file = joinpath(path, "temp.csm")
        open(temp_file, "w") do outf
            for line in eachline(g.csm_file)
                stripped = strip(line)
                if startswith(stripped, "despmtr")
                    parts = split(stripped)
                    pname = parts[2]
                    @assert pname in p.design_params.names
                    index = findfirst(pname, p.design_params.names)
                    # parts[3] = @sprintf("%.16g", param_vals[counter])
                    parts[3] = @sprintf("%.16g", param_vals[index])
                    counter = counter + 1
                    println(outf, join(parts, " "))
                else
                    println(outf, line)
                end
            end
        end
    end
end

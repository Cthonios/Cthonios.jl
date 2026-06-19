struct var"union (unnamed at EngSketchPad/include/egadsTypes.h:133:3)"
    data::NTuple{8, UInt8}
end

function Base.getproperty(x::Ptr{var"union (unnamed at EngSketchPad/include/egadsTypes.h:133:3)"}, f::Symbol)
    f === :integer && return Ptr{Cint}(x + 0)
    f === :integers && return Ptr{Ptr{Cint}}(x + 0)
    f === :real && return Ptr{Cdouble}(x + 0)
    f === :reals && return Ptr{Ptr{Cdouble}}(x + 0)
    f === :string && return Ptr{Ptr{Cchar}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"union (unnamed at EngSketchPad/include/egadsTypes.h:133:3)", f::Symbol)
    r = Ref{var"union (unnamed at EngSketchPad/include/egadsTypes.h:133:3)"}(x)
    ptr = Base.unsafe_convert(Ptr{var"union (unnamed at EngSketchPad/include/egadsTypes.h:133:3)"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"union (unnamed at EngSketchPad/include/egadsTypes.h:133:3)"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct egAttr
    data::NTuple{24, UInt8}
end

function Base.getproperty(x::Ptr{egAttr}, f::Symbol)
    f === :name && return Ptr{Ptr{Cchar}}(x + 0)
    f === :type && return Ptr{Cint}(x + 8)
    f === :length && return Ptr{Cint}(x + 12)
    f === :vals && return Ptr{var"union (unnamed at EngSketchPad/include/egadsTypes.h:133:3)"}(x + 16)
    return getfield(x, f)
end

function Base.getproperty(x::egAttr, f::Symbol)
    r = Ref{egAttr}(x)
    ptr = Base.unsafe_convert(Ptr{egAttr}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{egAttr}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct egAttrSeq
    root::Ptr{Cchar}
    nSeq::Cint
    attrSeq::Ptr{Cint}
end

struct egAttrs
    nattrs::Cint
    attrs::Ptr{egAttr}
    nseqs::Cint
    seqs::Ptr{egAttrSeq}
end

struct egObject
    magicnumber::Cint
    oclass::Cshort
    mtype::Cshort
    attrs::Ptr{Cvoid}
    blind::Ptr{Cvoid}
    topObj::Ptr{egObject}
    tref::Ptr{egObject}
    prev::Ptr{egObject}
    next::Ptr{egObject}
end

const ego = Ptr{egObject}

struct egCntxt
    outLevel::Cint
    fixedKnots::Cint
    fullAttrs::Cint
    tess::NTuple{2, Cdouble}
    signature::Ptr{Ptr{Cchar}}
    usrPtr::Ptr{Cvoid}
    threadID::Clong
    mutex::Ptr{Cvoid}
    pool::Ptr{egObject}
    last::Ptr{egObject}
end

struct egFconn
    index::Cint
    nface::Cint
    faces::Ptr{Cint}
    tric::Ptr{Cint}
end

struct egBary
    tri::Cint
    w::NTuple{2, Cdouble}
end

struct egTess1D
    obj::Ptr{egObject}
    nodes::NTuple{2, Cint}
    faces::NTuple{2, egFconn}
    xyz::Ptr{Cdouble}
    t::Ptr{Cdouble}
    _global::Ptr{Cint}
    npts::Cint
end

struct egPatch
    ipts::Ptr{Cint}
    bounds::Ptr{Cint}
    nu::Cint
    nv::Cint
end

struct egTess2D
    mKnots::Ptr{egObject}
    xyz::Ptr{Cdouble}
    uv::Ptr{Cdouble}
    _global::Ptr{Cint}
    ptype::Ptr{Cint}
    pindex::Ptr{Cint}
    bary::Ptr{egBary}
    frame::Ptr{Cint}
    frlps::Ptr{Cint}
    tris::Ptr{Cint}
    tric::Ptr{Cint}
    patch::Ptr{egPatch}
    npts::Cint
    nframe::Cint
    nfrlps::Cint
    ntris::Cint
    npatch::Cint
    tfi::Cint
end

struct egTessel
    src::Ptr{egObject}
    xyzs::Ptr{Cdouble}
    tess1d::Ptr{egTess1D}
    tess2d::Ptr{egTess2D}
    globals::Ptr{Cint}
    params::NTuple{6, Cdouble}
    tparam::NTuple{2, Cdouble}
    nGlobal::Cint
    nEdge::Cint
    nFace::Cint
    nu::Cint
    nv::Cint
    done::Cint
end

struct egEMap
    nobjs::Cint
    objs::Ptr{Ptr{egObject}}
end

struct egEdVert
    edge::Ptr{egObject}
    curve::Cint
    npts::Cint
    ts::Ptr{Cdouble}
    dstart::NTuple{3, Cdouble}
    dend::NTuple{3, Cdouble}
end

struct egEEseg
    iedge::Cint
    sense::Cint
    nstart::Ptr{egObject}
    tstart::Cdouble
    tend::Cdouble
end

struct egEEdge
    sedges::Ptr{egEdVert}
    nsegs::Cint
    segs::Ptr{egEEseg}
    trange::NTuple{2, Cdouble}
    nodes::NTuple{2, Ptr{egObject}}
end

struct egEdgeUV
    edge::Ptr{egObject}
    sense::Cint
    npts::Cint
    iuv::Ptr{Cint}
end

struct egELoop
    eedges::egEMap
    senses::Ptr{Cint}
    area::Cdouble
    nedge::Cint
    edgeUVs::Ptr{egEdgeUV}
end

struct egEPatch
    face::Ptr{egObject}
    tol::Cdouble
    start::Cint
    nuvs::Cint
    ndeflect::Cint
    ntris::Cint
    uvtris::Ptr{Cint}
    uvtric::Ptr{Cint}
    uvs::Ptr{Cdouble}
    deflect::Ptr{Cdouble}
end

struct egEFace
    sedges::Ptr{egEdVert}
    npatch::Cint
    patches::Ptr{egEPatch}
    eloops::egEMap
    senses::Ptr{Cint}
    trmap::Ptr{Cint}
    uvmap::Ptr{Cvoid}
    range::NTuple{4, Cdouble}
    last::Cint
end

struct egEShell
    efaces::egEMap
end

struct egEBody
    ref::Ptr{egObject}
    eedges::egEMap
    eloops::egEMap
    efaces::egEMap
    eshells::egEMap
    senses::Ptr{Cint}
    angle::Cdouble
    done::Cint
    nedge::Cint
    edges::Ptr{egEdVert}
end

function caps_revision(major, minor)
    ccall((:caps_revision, "libcaps"), Cvoid, (Ptr{Cint}, Ptr{Cint}), major, minor)
end

struct capsOwn
    index::Cint
    pname::Ptr{Cchar}
    pID::Ptr{Cchar}
    user::Ptr{Cchar}
    datetime::NTuple{6, Cshort}
    sNum::Culong
end

struct capsObject
    magicnumber::Cint
    type::Cint
    subtype::Cint
    delMark::Cint
    name::Ptr{Cchar}
    attrs::Ptr{egAttrs}
    blind::Ptr{Cvoid}
    flist::Ptr{Cvoid}
    nHistory::Cint
    history::Ptr{capsOwn}
    last::capsOwn
    parent::Ptr{capsObject}
end

const capsObj = Ptr{capsObject}

@cenum capsoType::Int32 begin
    BODIES = -2
    ATTRIBUTES = -1
    UNUSED = 0
    PROBLEM = 1
    VALUE = 2
    ANALYSIS = 3
    BOUND = 4
    VERTEXSET = 5
    DATASET = 6
end

@cenum capssType::UInt32 begin
    NONE = 0
    STATIC = 1
    PARAMETRIC = 2
    GEOMETRYIN = 3
    GEOMETRYOUT = 4
    PARAMETER = 5
    USER = 6
    ANALYSISIN = 7
    ANALYSISOUT = 8
    CONNECTED = 9
    UNCONNECTED = 10
    ANALYSISDYNO = 11
end

function caps_info(object, name, type, subtype, link, parent, last)
    ccall((:caps_info, "libcaps"), Cint, (capsObj, Ptr{Ptr{Cchar}}, Ptr{capsoType}, Ptr{capssType}, Ptr{capsObj}, Ptr{capsObj}, Ptr{capsOwn}), object, name, type, subtype, link, parent, last)
end

struct capsError
    errObj::Ptr{capsObject}
    eType::Cint
    index::Cint
    nLines::Cint
    lines::Ptr{Ptr{Cchar}}
end

struct capsErrs
    nError::Cint
    errors::Ptr{capsError}
end

function caps_size(object, type, stype, size, nErr, errors)
    ccall((:caps_size, "libcaps"), Cint, (capsObj, capsoType, capssType, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), object, type, stype, size, nErr, errors)
end

function caps_childByIndex(object, type, stype, index, child)
    ccall((:caps_childByIndex, "libcaps"), Cint, (capsObj, capsoType, capssType, Cint, Ptr{capsObj}), object, type, stype, index, child)
end

function caps_childByName(object, typ, styp, name, child, nErr, errors)
    ccall((:caps_childByName, "libcaps"), Cint, (capsObj, capsoType, capssType, Ptr{Cchar}, Ptr{capsObj}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), object, typ, styp, name, child, nErr, errors)
end

function caps_bodyByIndex(pobject, index, body, units)
    ccall((:caps_bodyByIndex, "libcaps"), Cint, (capsObj, Cint, Ptr{ego}, Ptr{Ptr{Cchar}}), pobject, index, body, units)
end

function caps_ownerInfo(pobject, owner, phase, pname, pID, userID, nLines, lines, datetime, sNum)
    ccall((:caps_ownerInfo, "libcaps"), Cint, (capsObj, capsOwn, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cshort}, Ptr{Culong}), pobject, owner, phase, pname, pID, userID, nLines, lines, datetime, sNum)
end

function caps_getHistory(object, nHist, history)
    ccall((:caps_getHistory, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Ptr{capsOwn}}), object, nHist, history)
end

function caps_markForDelete(object)
    ccall((:caps_markForDelete, "libcaps"), Cint, (capsObj,), object)
end

function caps_errorInfo(errs, eIndex, errObj, eType, nLines, lines)
    ccall((:caps_errorInfo, "libcaps"), Cint, (Ptr{capsErrs}, Cint, Ptr{capsObj}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}), errs, eIndex, errObj, eType, nLines, lines)
end

function caps_freeError(errs)
    ccall((:caps_freeError, "libcaps"), Cint, (Ptr{capsErrs},), errs)
end

function caps_printErrors(fp, nErr, errs)
    ccall((:caps_printErrors, "libcaps"), Cint, (Ptr{Cint}, Cint, Ptr{capsErrs}), fp, nErr, errs)
end

struct var"union (unnamed at EngSketchPad/include/capsTypes.h:291:3)"
    data::NTuple{8, UInt8}
end

function Base.getproperty(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:291:3)"}, f::Symbol)
    f === :integer && return Ptr{Cint}(x + 0)
    f === :integers && return Ptr{Ptr{Cint}}(x + 0)
    f === :real && return Ptr{Cdouble}(x + 0)
    f === :reals && return Ptr{Ptr{Cdouble}}(x + 0)
    f === :string && return Ptr{Ptr{Cchar}}(x + 0)
    f === :tuple && return Ptr{Ptr{capsTuple}}(x + 0)
    f === :AIMptr && return Ptr{Ptr{Cvoid}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"union (unnamed at EngSketchPad/include/capsTypes.h:291:3)", f::Symbol)
    r = Ref{var"union (unnamed at EngSketchPad/include/capsTypes.h:291:3)"}(x)
    ptr = Base.unsafe_convert(Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:291:3)"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:291:3)"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct var"union (unnamed at EngSketchPad/include/capsTypes.h:300:3)"
    data::NTuple{16, UInt8}
end

function Base.getproperty(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:300:3)"}, f::Symbol)
    f === :ilims && return Ptr{NTuple{2, Cint}}(x + 0)
    f === :dlims && return Ptr{NTuple{2, Cdouble}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"union (unnamed at EngSketchPad/include/capsTypes.h:300:3)", f::Symbol)
    r = Ref{var"union (unnamed at EngSketchPad/include/capsTypes.h:300:3)"}(x)
    ptr = Base.unsafe_convert(Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:300:3)"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:300:3)"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct capsDeriv
    name::Ptr{Cchar}
    len_wrt::Cint
    deriv::Ptr{Cdouble}
end

struct capsValue
    data::NTuple{144, UInt8}
end

function Base.getproperty(x::Ptr{capsValue}, f::Symbol)
    f === :type && return Ptr{Cint}(x + 0)
    f === :length && return Ptr{Cint}(x + 4)
    f === :dim && return Ptr{Cint}(x + 8)
    f === :nrow && return Ptr{Cint}(x + 12)
    f === :ncol && return Ptr{Cint}(x + 16)
    f === :lfixed && return Ptr{Cint}(x + 20)
    f === :sfixed && return Ptr{Cint}(x + 24)
    f === :nullVal && return Ptr{Cint}(x + 28)
    f === :index && return Ptr{Cint}(x + 32)
    f === :pIndex && return Ptr{Cint}(x + 36)
    f === :gInType && return Ptr{Cint}(x + 40)
    f === :vals && return Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:291:3)"}(x + 48)
    f === :limits && return Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:300:3)"}(x + 56)
    f === :lims && return Ptr{Ptr{Cvoid}}(x + 72)
    f === :units && return Ptr{Ptr{Cchar}}(x + 80)
    f === :meshWriter && return Ptr{Ptr{Cchar}}(x + 88)
    f === :link && return Ptr{Ptr{capsObject}}(x + 96)
    f === :linkMethod && return Ptr{Cint}(x + 104)
    f === :partial && return Ptr{Ptr{Cint}}(x + 112)
    f === :nderiv && return Ptr{Cint}(x + 120)
    f === :derivs && return Ptr{Ptr{capsDeriv}}(x + 128)
    f === :stepSize && return Ptr{Ptr{Cdouble}}(x + 136)
    return getfield(x, f)
end

function Base.getproperty(x::capsValue, f::Symbol)
    r = Ref{capsValue}(x)
    ptr = Base.unsafe_convert(Ptr{capsValue}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{capsValue}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function caps_freeValue(value)
    ccall((:caps_freeValue, "libcaps"), Cvoid, (Ptr{capsValue},), value)
end

function caps_writeParameters(pobject, filename)
    ccall((:caps_writeParameters, "libcaps"), Cint, (capsObj, Ptr{Cchar}), pobject, filename)
end

function caps_readParameters(pobject, filename)
    ccall((:caps_readParameters, "libcaps"), Cint, (capsObj, Ptr{Cchar}), pobject, filename)
end

function caps_writeGeometry(object, flag, filename, nErr, errors)
    ccall((:caps_writeGeometry, "libcaps"), Cint, (capsObj, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), object, flag, filename, nErr, errors)
end

function caps_attrByName(cobj, name, attr)
    ccall((:caps_attrByName, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Ptr{capsObj}), cobj, name, attr)
end

function caps_attrByIndex(cobj, index, attr)
    ccall((:caps_attrByIndex, "libcaps"), Cint, (capsObj, Cint, Ptr{capsObj}), cobj, index, attr)
end

function caps_setAttr(cobj, name, aval)
    ccall((:caps_setAttr, "libcaps"), Cint, (capsObj, Ptr{Cchar}, capsObj), cobj, name, aval)
end

function caps_deleteAttr(cobj, name)
    ccall((:caps_deleteAttr, "libcaps"), Cint, (capsObj, Ptr{Cchar}), cobj, name)
end

function caps_phaseState(prName, phName, bitFlag)
    ccall((:caps_phaseState, "libcaps"), Cint, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}), prName, phName, bitFlag)
end

function caps_phaseNewCSM(prName, phName, csm)
    ccall((:caps_phaseNewCSM, "libcaps"), Cint, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}), prName, phName, csm)
end

function caps_journalState(pobject)
    ccall((:caps_journalState, "libcaps"), Cint, (capsObj,), pobject)
end

function caps_open(prName, phName, flag, ptr, outLevel, pobject, nErr, errors)
    ccall((:caps_open, "libcaps"), Cint, (Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Cvoid}, Cint, Ptr{capsObj}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), prName, phName, flag, ptr, outLevel, pobject, nErr, errors)
end

function caps_brokenLink(callBack)
    ccall((:caps_brokenLink, "libcaps"), Cint, (Ptr{Cvoid},), callBack)
end

function caps_close(pobject, complete, phName)
    ccall((:caps_close, "libcaps"), Cint, (capsObj, Cint, Ptr{Cchar}), pobject, complete, phName)
end

function caps_outLevel(pobject, outLevel)
    ccall((:caps_outLevel, "libcaps"), Cint, (capsObj, Cint), pobject, outLevel)
end

function caps_getRootPath(pobject, fullPath)
    ccall((:caps_getRootPath, "libcaps"), Cint, (capsObj, Ptr{Ptr{Cchar}}), pobject, fullPath)
end

function caps_intentPhrase(pobject, nLines, lines)
    ccall((:caps_intentPhrase, "libcaps"), Cint, (capsObj, Cint, Ptr{Ptr{Cchar}}), pobject, nLines, lines)
end

function caps_debug(pobject)
    ccall((:caps_debug, "libcaps"), Cint, (capsObj,), pobject)
end

function caps_queryAnalysis(pobj, aname, nIn, nOut, execution)
    ccall((:caps_queryAnalysis, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), pobj, aname, nIn, nOut, execution)
end

function caps_getBodies(aobject, nBody, bodies, nErr, errors)
    ccall((:caps_getBodies, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), aobject, nBody, bodies, nErr, errors)
end

function caps_execute(object, state, nErr, errors)
    ccall((:caps_execute, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), object, state, nErr, errors)
end

function caps_getInput(pobj, aname, index, ainame, defaults)
    ccall((:caps_getInput, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Cint, Ptr{Ptr{Cchar}}, Ptr{capsValue}), pobj, aname, index, ainame, defaults)
end

function caps_getOutput(pobj, aname, index, aoname, form)
    ccall((:caps_getOutput, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Cint, Ptr{Ptr{Cchar}}, Ptr{capsValue}), pobj, aname, index, aoname, form)
end

function caps_AIMbackdoor(aobject, JSONin, JSONout)
    ccall((:caps_AIMbackdoor, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Ptr{Ptr{Cchar}}), aobject, JSONin, JSONout)
end

function caps_makeAnalysis(pobj, anam, nam, unitSys, intents, exec, aobject, nErr, errors)
    ccall((:caps_makeAnalysis, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{capsObj}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), pobj, anam, nam, unitSys, intents, exec, aobject, nErr, errors)
end

function caps_dupAnalysis(from, name, aobj)
    ccall((:caps_dupAnalysis, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Ptr{capsObj}), from, name, aobj)
end

function caps_dirtyAnalysis(pobj, nAobj, aobjs)
    ccall((:caps_dirtyAnalysis, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Ptr{capsObj}}), pobj, nAobj, aobjs)
end

function caps_analysisInfo(aobject, apath, unitSys, major, minor, intents, nField, fnames, ranks, fInOut, execute, status)
    ccall((:caps_analysisInfo, "libcaps"), Cint, (capsObj, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Cint}), aobject, apath, unitSys, major, minor, intents, nField, fnames, ranks, fInOut, execute, status)
end

function caps_preAnalysis(aobject, nErr, errors)
    ccall((:caps_preAnalysis, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Ptr{capsErrs}}), aobject, nErr, errors)
end

function caps_system(aobject, rpath, command)
    ccall((:caps_system, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Ptr{Cchar}), aobject, rpath, command)
end

function caps_postAnalysis(aobject, nErr, errors)
    ccall((:caps_postAnalysis, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Ptr{capsErrs}}), aobject, nErr, errors)
end

function caps_getTessels(aobject, nTessel, tessels, nErr, errors)
    ccall((:caps_getTessels, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), aobject, nTessel, tessels, nErr, errors)
end

function caps_makeBound(pobject, dim, bname, bobj)
    ccall((:caps_makeBound, "libcaps"), Cint, (capsObj, Cint, Ptr{Cchar}, Ptr{capsObj}), pobject, dim, bname, bobj)
end

@cenum capsState::Int32 begin
    MultipleError = -2
    Open = -1
    Empty = 0
    Single = 1
    Multiple = 2
end

function caps_boundInfo(object, state, dim, plims)
    ccall((:caps_boundInfo, "libcaps"), Cint, (capsObj, Ptr{capsState}, Ptr{Cint}, Ptr{Cdouble}), object, state, dim, plims)
end

function caps_closeBound(bobject)
    ccall((:caps_closeBound, "libcaps"), Cint, (capsObj,), bobject)
end

function caps_makeVertexSet(bobject, aobject, vname, vobj, nErr, errors)
    ccall((:caps_makeVertexSet, "libcaps"), Cint, (capsObj, capsObj, Ptr{Cchar}, Ptr{capsObj}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), bobject, aobject, vname, vobj, nErr, errors)
end

function caps_vertexSetInfo(vobject, nGpts, nDpts, bobj, aobj)
    ccall((:caps_vertexSetInfo, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Cint}, Ptr{capsObj}, Ptr{capsObj}), vobject, nGpts, nDpts, bobj, aobj)
end

function caps_outputVertexSet(vobject, filename)
    ccall((:caps_outputVertexSet, "libcaps"), Cint, (capsObj, Ptr{Cchar}), vobject, filename)
end

function caps_fillUnVertexSet(vobject, npts, xyzs)
    ccall((:caps_fillUnVertexSet, "libcaps"), Cint, (capsObj, Cint, Ptr{Cdouble}), vobject, npts, xyzs)
end

@cenum capsfType::UInt32 begin
    FieldIn = 0
    FieldOut = 1
    GeomSens = 2
    TessSens = 3
    User = 4
    BuiltIn = 5
end

function caps_makeDataSet(vobject, dname, ftype, rank, dobj, nErr, errors)
    ccall((:caps_makeDataSet, "libcaps"), Cint, (capsObj, Ptr{Cchar}, capsfType, Cint, Ptr{capsObj}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), vobject, dname, ftype, rank, dobj, nErr, errors)
end

@cenum capsdMethod::UInt32 begin
    Interpolate = 0
    Conserve = 1
end

function caps_dataSetInfo(dobject, ftype, link, dmeth)
    ccall((:caps_dataSetInfo, "libcaps"), Cint, (capsObj, Ptr{capsfType}, Ptr{capsObj}, Ptr{capsdMethod}), dobject, ftype, link, dmeth)
end

function caps_linkDataSet(link, method, target, nErr, errors)
    ccall((:caps_linkDataSet, "libcaps"), Cint, (capsObj, capsdMethod, capsObj, Ptr{Cint}, Ptr{Ptr{capsErrs}}), link, method, target, nErr, errors)
end

function caps_initDataSet(dobject, rank, data, nErr, errors)
    ccall((:caps_initDataSet, "libcaps"), Cint, (capsObj, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), dobject, rank, data, nErr, errors)
end

function caps_setData(dobject, nverts, rank, data, units, nErr, errors)
    ccall((:caps_setData, "libcaps"), Cint, (capsObj, Cint, Cint, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), dobject, nverts, rank, data, units, nErr, errors)
end

function caps_getData(dobject, npts, rank, data, units, nErr, errors)
    ccall((:caps_getData, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), dobject, npts, rank, data, units, nErr, errors)
end

function caps_getDataSets(bobject, dname, nobj, dobjs)
    ccall((:caps_getDataSets, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Ptr{Cint}, Ptr{Ptr{capsObj}}), bobject, dname, nobj, dobjs)
end

function caps_getTriangles(vobject, nGtris, gtris, nGsegs, gsegs, nDtris, dtris, nDsegs, dsegs)
    ccall((:caps_getTriangles, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}), vobject, nGtris, gtris, nGsegs, gsegs, nDtris, dtris, nDsegs, dsegs)
end

@cenum capsvType::Int32 begin
    Doubles = -2
    Integers = -1
    Boolean = 0
    Integer = 1
    Double = 2
    mString = 3
    Tuple = 4
    Pointer = 5
    DoubleDeriv = 6
    PointerMesh = 7
end

function caps_getValue(object, vtype, nrow, ncol, data, partial, units, nErr, errors)
    ccall((:caps_getValue, "libcaps"), Cint, (capsObj, Ptr{capsvType}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), object, vtype, nrow, ncol, data, partial, units, nErr, errors)
end

function caps_makeValue(pobject, vname, stype, vtype, nrow, ncol, data, partial, units, vobj)
    ccall((:caps_makeValue, "libcaps"), Cint, (capsObj, Ptr{Cchar}, capssType, capsvType, Cint, Cint, Ptr{Cvoid}, Ptr{Cint}, Ptr{Cchar}, Ptr{capsObj}), pobject, vname, stype, vtype, nrow, ncol, data, partial, units, vobj)
end

function caps_setValue(object, vtype, nrow, ncol, data, partial, units, nErr, errors)
    ccall((:caps_setValue, "libcaps"), Cint, (capsObj, capsvType, Cint, Cint, Ptr{Cvoid}, Ptr{Cint}, Ptr{Cchar}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), object, vtype, nrow, ncol, data, partial, units, nErr, errors)
end

function caps_getLimits(object, vtype, limits, units)
    ccall((:caps_getLimits, "libcaps"), Cint, (capsObj, Ptr{capsvType}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cchar}}), object, vtype, limits, units)
end

function caps_setLimits(object, vtype, limits, units, nErr, errors)
    ccall((:caps_setLimits, "libcaps"), Cint, (capsObj, capsvType, Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), object, vtype, limits, units, nErr, errors)
end

function caps_getValueSize(object, nrow, ncol)
    ccall((:caps_getValueSize, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Cint}), object, nrow, ncol)
end

@cenum capsFixed::UInt32 begin
    Change = 0
    Fixed = 1
end

@cenum capsNull::UInt32 begin
    NotAllowed = 0
    NotNull = 1
    IsNull = 2
    IsPartial = 3
end

function caps_getValueProps(object, dim, gInType, lfix, sfix, ntype)
    ccall((:caps_getValueProps, "libcaps"), Cint, (capsObj, Ptr{Cint}, Ptr{Cint}, Ptr{capsFixed}, Ptr{capsFixed}, Ptr{capsNull}), object, dim, gInType, lfix, sfix, ntype)
end

function caps_setValueProps(object, dim, lfixed, sfixed, ntype, nErr, errors)
    ccall((:caps_setValueProps, "libcaps"), Cint, (capsObj, Cint, capsFixed, capsFixed, capsNull, Ptr{Cint}, Ptr{Ptr{capsErrs}}), object, dim, lfixed, sfixed, ntype, nErr, errors)
end

function caps_convertValue(value, inVal, inUnit, outVal)
    ccall((:caps_convertValue, "libcaps"), Cint, (capsObj, Cdouble, Ptr{Cchar}, Ptr{Cdouble}), value, inVal, inUnit, outVal)
end

@cenum capstMethod::UInt32 begin
    Copy = 0
    Integrate = 1
    Average = 2
end

function caps_transferValues(source, method, target, nErr, errors)
    ccall((:caps_transferValues, "libcaps"), Cint, (capsObj, capstMethod, capsObj, Ptr{Cint}, Ptr{Ptr{capsErrs}}), source, method, target, nErr, errors)
end

function caps_linkValue(link, method, target, nErr, errors)
    ccall((:caps_linkValue, "libcaps"), Cint, (capsObj, capstMethod, capsObj, Ptr{Cint}, Ptr{Ptr{capsErrs}}), link, method, target, nErr, errors)
end

function caps_hasDeriv(vobj, nderiv, names, nErr, errors)
    ccall((:caps_hasDeriv, "libcaps"), Cint, (Ptr{capsObject}, Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), vobj, nderiv, names, nErr, errors)
end

function caps_getDeriv(value, name, len, len_wrt, deriv, nErr, errors)
    ccall((:caps_getDeriv, "libcaps"), Cint, (capsObj, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Cint}, Ptr{Ptr{capsErrs}}), value, name, len, len_wrt, deriv, nErr, errors)
end

function caps_setStepSize(object, steps)
    ccall((:caps_setStepSize, "libcaps"), Cint, (capsObj, Ptr{Cdouble}), object, steps)
end

function caps_getStepSize(object, steps)
    ccall((:caps_getStepSize, "libcaps"), Cint, (capsObj, Ptr{Ptr{Cdouble}}), object, steps)
end

function caps_convert(count, inUnit, inVal, outUnit, outVal)
    ccall((:caps_convert, "libcaps"), Cint, (Cint, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cdouble}), count, inUnit, inVal, outUnit, outVal)
end

function caps_unitParse(unit)
    ccall((:caps_unitParse, "libcaps"), Cint, (Ptr{Cchar},), unit)
end

function caps_unitConvertible(unit1, unit2)
    ccall((:caps_unitConvertible, "libcaps"), Cint, (Ptr{Cchar}, Ptr{Cchar}), unit1, unit2)
end

function caps_unitCompare(unit1, unit2, compare)
    ccall((:caps_unitCompare, "libcaps"), Cint, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}), unit1, unit2, compare)
end

function caps_unitMultiply(inUnits1, inUnits2, outUnits)
    ccall((:caps_unitMultiply, "libcaps"), Cint, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), inUnits1, inUnits2, outUnits)
end

function caps_unitDivide(inUnits1, inUnits2, outUnits)
    ccall((:caps_unitDivide, "libcaps"), Cint, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), inUnits1, inUnits2, outUnits)
end

function caps_unitInvert(inUnit, outUnits)
    ccall((:caps_unitInvert, "libcaps"), Cint, (Ptr{Cchar}, Ptr{Ptr{Cchar}}), inUnit, outUnits)
end

function caps_unitRaise(inUnit, power, outUnits)
    ccall((:caps_unitRaise, "libcaps"), Cint, (Ptr{Cchar}, Cint, Ptr{Ptr{Cchar}}), inUnit, power, outUnits)
end

function caps_unitOffset(inUnit, offset, outUnits)
    ccall((:caps_unitOffset, "libcaps"), Cint, (Ptr{Cchar}, Cdouble, Ptr{Ptr{Cchar}}), inUnit, offset, outUnits)
end

# no prototype is found for this function at caps.h:375:3, please use with caution
function caps_externSignal()
    ccall((:caps_externSignal, "libcaps"), Cvoid, ())
end

# no prototype is found for this function at caps.h:378:3, please use with caution
function caps_rmLock()
    ccall((:caps_rmLock, "libcaps"), Cvoid, ())
end

function caps_printObjects(pobject, object, indent)
    ccall((:caps_printObjects, "libcaps"), Cvoid, (capsObj, capsObj, Cint), pobject, object, indent)
end

@cenum capsoFlag::UInt32 begin
    oFileName = 0
    oMODL = 1
    oEGO = 2
    oPhaseName = 3
    oContinue = 4
    oPNewCSM = 5
    oPNnoDel = 6
    oReadOnly = 7
end

@cenum capseType::Int32 begin
    CONTINUATION = -1
    CINFO = 0
    CWARN = 1
    CERROR = 2
    CSTAT = 3
end

@cenum capsjType::UInt32 begin
    jInteger = 0
    jDouble = 1
    jString = 2
    jStrings = 3
    jTuple = 4
    jPointer = 5
    jPtrFree = 6
    jObject = 7
    jObjs = 8
    jErr = 9
    jOwn = 10
    jOwns = 11
    jEgos = 12
    jFile = 13
end

@cenum capsBoolean::UInt32 begin
    False = 0
    True = 1
end

@cenum capsvDim::UInt32 begin
    Scalar = 0
    mVector = 1
    Array2D = 2
end

struct capsTuple
    name::Ptr{Cchar}
    value::Ptr{Cchar}
end

struct capsAprx1D
    nrank::Cint
    periodic::Cint
    nts::Cint
    interp::Ptr{Cdouble}
    trange::NTuple{2, Cdouble}
    ntm::Cint
    tmap::Ptr{Cdouble}
end

struct capsAprx2D
    nrank::Cint
    periodic::Cint
    nus::Cint
    nvs::Cint
    interp::Ptr{Cdouble}
    urange::NTuple{2, Cdouble}
    vrange::NTuple{2, Cdouble}
    num::Cint
    nvm::Cint
    uvmap::Ptr{Cdouble}
end

struct capsEleType
    nref::Cint
    ndata::Cint
    nmat::Cint
    ntri::Cint
    gst::Ptr{Cdouble}
    dst::Ptr{Cdouble}
    matst::Ptr{Cdouble}
    tris::Ptr{Cint}
    nseg::Cint
    segs::Ptr{Cint}
end

struct var"union (unnamed at EngSketchPad/include/capsTypes.h:164:3)"
    data::NTuple{8, UInt8}
end

function Base.getproperty(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:164:3)"}, f::Symbol)
    f === :tq && return Ptr{NTuple{2, Cint}}(x + 0)
    f === :poly && return Ptr{Ptr{Cint}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"union (unnamed at EngSketchPad/include/capsTypes.h:164:3)", f::Symbol)
    r = Ref{var"union (unnamed at EngSketchPad/include/capsTypes.h:164:3)"}(x)
    ptr = Base.unsafe_convert(Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:164:3)"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:164:3)"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct capsElement
    data::NTuple{32, UInt8}
end

function Base.getproperty(x::Ptr{capsElement}, f::Symbol)
    f === :tIndex && return Ptr{Cint}(x + 0)
    f === :eIndex && return Ptr{Cint}(x + 4)
    f === :gIndices && return Ptr{Ptr{Cint}}(x + 8)
    f === :dIndices && return Ptr{Ptr{Cint}}(x + 16)
    f === :eTris && return Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:164:3)"}(x + 24)
    return getfield(x, f)
end

function Base.getproperty(x::capsElement, f::Symbol)
    r = Ref{capsElement}(x)
    ptr = Base.unsafe_convert(Ptr{capsElement}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{capsElement}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct capsBodyDiscr
    tess::ego
    nElems::Cint
    elems::Ptr{capsElement}
    gIndices::Ptr{Cint}
    dIndices::Ptr{Cint}
    poly::Ptr{Cint}
    globalOffset::Cint
end

struct capsDiscr
    dim::Cint
    instStore::Ptr{Cvoid}
    aInfo::Ptr{Cvoid}
    nVerts::Cint
    verts::Ptr{Cdouble}
    celem::Ptr{Cint}
    nDtris::Cint
    dtris::Ptr{Cint}
    nDsegs::Cint
    dsegs::Ptr{Cint}
    nPoints::Cint
    nTypes::Cint
    types::Ptr{capsEleType}
    nBodys::Cint
    bodys::Ptr{capsBodyDiscr}
    tessGlobal::Ptr{Cint}
    ptrm::Ptr{Cvoid}
end

# struct capsJrnl
#     type::Cint
#     num::Cint
#     length::Cint
#     members::var"union (unnamed at EngSketchPad/include/capsTypes.h:323:3)"
# end

struct var"union (unnamed at EngSketchPad/include/capsTypes.h:346:3)"
    data::NTuple{56, UInt8}
end

function Base.getproperty(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:346:3)"}, f::Symbol)
    f === :tuple && return Ptr{Ptr{capsTuple}}(x + 0)
    f === :strings && return Ptr{Ptr{Ptr{Cchar}}}(x + 0)
    f === :pointer && return Ptr{Ptr{Cvoid}}(x + 0)
    f === :own && return Ptr{capsOwn}(x + 0)
    f === :owns && return Ptr{Ptr{capsOwn}}(x + 0)
    f === :model && return Ptr{ego}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"union (unnamed at EngSketchPad/include/capsTypes.h:346:3)", f::Symbol)
    r = Ref{var"union (unnamed at EngSketchPad/include/capsTypes.h:346:3)"}(x)
    ptr = Base.unsafe_convert(Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:346:3)"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:346:3)"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct capsFList
    data::NTuple{80, UInt8}
end

function Base.getproperty(x::Ptr{capsFList}, f::Symbol)
    f === :type && return Ptr{Cint}(x + 0)
    f === :num && return Ptr{Cint}(x + 4)
    f === :member && return Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:346:3)"}(x + 8)
    f === :sNum && return Ptr{Culong}(x + 64)
    f === :next && return Ptr{Ptr{capsFList}}(x + 72)
    return getfield(x, f)
end

function Base.getproperty(x::capsFList, f::Symbol)
    r = Ref{capsFList}(x)
    ptr = Base.unsafe_convert(Ptr{capsFList}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{capsFList}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

# typedef int ( * aimI ) ( int , /*@null@*/ const char * , /*@null@*/ void * , /*@null@*/ void * * , int * , int * , int * , int * , int * , char * * * , int * * , int * * )
const aimI = Ptr{Cvoid}

# typedef int ( * aimD ) ( char * , capsDiscr * )
const aimD = Ptr{Cvoid}

# typedef void ( * aimF ) ( void * )
const aimF = Ptr{Cvoid}

# typedef int ( * aimL ) ( capsDiscr * , double * , double * , int * , int * , double * )
const aimL = Ptr{Cvoid}

# typedef int ( * aimIn ) ( /*@null@*/ void * , /*@null@*/ void * , int , char * * , capsValue * )
const aimIn = Ptr{Cvoid}

# typedef int ( * aimU ) ( /*@null@*/ void * , void * , /*@null@*/ capsValue * )
const aimU = Ptr{Cvoid}

# typedef int ( * aimA ) ( /*@null@*/ const void * , void * , /*@null@*/ capsValue * )
const aimA = Ptr{Cvoid}

# typedef int ( * aimEx ) ( /*@null@*/ const void * , void * , int * )
const aimEx = Ptr{Cvoid}

# typedef int ( * aimPo ) ( /*@null@*/ void * , void * , int , /*@null@*/ capsValue * )
const aimPo = Ptr{Cvoid}

# typedef int ( * aimO ) ( /*@null@*/ void * , /*@null@*/ void * , int , char * * , capsValue * )
const aimO = Ptr{Cvoid}

# typedef int ( * aimC ) ( /*@null@*/ void * , void * , int , capsValue * )
const aimC = Ptr{Cvoid}

# typedef int ( * aimT ) ( capsDiscr * , const char * , int , int , double * , char * * )
const aimT = Ptr{Cvoid}

# typedef int ( * aimP ) ( capsDiscr * , const char * , int , int , double * , int , double * , double * )
const aimP = Ptr{Cvoid}

# typedef int ( * aimG ) ( capsDiscr * , const char * , int , int , int , /*@null@*/ double * , double * )
const aimG = Ptr{Cvoid}

# typedef int ( * aimDa ) ( /*@null@*/ void * , void * , const char * , enum capsvType * , int * , int * , int * , void * * , char * * )
const aimDa = Ptr{Cvoid}

# typedef int ( * aimBd ) ( /*@null@*/ void * , void * , const char * , char * * )
const aimBd = Ptr{Cvoid}

# typedef void ( * aimCU ) ( void * )
const aimCU = Ptr{Cvoid}

# typedef int ( * AIMwriter ) ( void * , void * )
const AIMwriter = Ptr{Cvoid}

# typedef const char * ( * AIMext ) ( void )
const AIMext = Ptr{Cvoid}

struct aimContext
    aim_nAnal::Cint
    aimName::NTuple{64, Ptr{Cchar}}
    aim_nInst::NTuple{64, Cint}
    aimDLL::NTuple{64, Ptr{Cvoid}}
    aimInit::NTuple{64, aimI}
    aimDiscr::NTuple{64, aimD}
    aimFreeD::NTuple{64, aimF}
    aimLoc::NTuple{64, aimL}
    aimInput::NTuple{64, aimIn}
    aimUState::NTuple{64, aimU}
    aimPAnal::NTuple{64, aimA}
    aimExec::NTuple{64, aimEx}
    aimPost::NTuple{64, aimPo}
    aimOutput::NTuple{64, aimO}
    aimCalc::NTuple{64, aimC}
    aimXfer::NTuple{64, aimT}
    aimIntrp::NTuple{64, aimP}
    aimIntrpBar::NTuple{64, aimP}
    aimIntgr::NTuple{64, aimG}
    aimIntgrBar::NTuple{64, aimG}
    aimBdoor::NTuple{64, aimBd}
    aimClean::NTuple{64, aimCU}
end

struct capsRegGIN
    name::Ptr{Cchar}
    index::Cint
    irow::Cint
    icol::Cint
end

struct capsPhrase
    phase::Ptr{Cchar}
    nLines::Cint
    lines::Ptr{Ptr{Cchar}}
end

struct capsProblem
    signature::Ptr{Ptr{Cchar}}
    mySelf::Ptr{capsObject}
    context::ego
    utsystem::Ptr{Cvoid}
    aimFPTR::aimContext
    root::Ptr{Cchar}
    phName::Ptr{Cchar}
    writer::capsOwn
    dbFlag::Cint
    stFlag::Cint
    jrnl::Ptr{Cint}
    outLevel::Cint
    funID::Cint
    modl::Ptr{Cvoid}
    DTime::Cdouble
    iPhrase::Cint
    nPhrase::Cint
    phrases::Ptr{capsPhrase}
    nParam::Cint
    params::Ptr{Ptr{capsObject}}
    nUser::Cint
    users::Ptr{Ptr{capsObject}}
    nGeomIn::Cint
    geomIn::Ptr{Ptr{capsObject}}
    nGeomOut::Cint
    geomOut::Ptr{Ptr{capsObject}}
    nAnalysis::Cint
    analysis::Ptr{Ptr{capsObject}}
    mBound::Cint
    nBound::Cint
    bounds::Ptr{Ptr{capsObject}}
    geometry::capsOwn
    nBodies::Cint
    bodies::Ptr{ego}
    lunits::Ptr{Ptr{Cchar}}
    nEGADSmdl::Cint
    nRegGIN::Cint
    regGIN::Ptr{capsRegGIN}
    nDesPmtr::Cint
    desPmtr::Ptr{Cint}
    sNum::Culong
    nFiles::Cint
    files::Ptr{Culong}
    jpos::Clong
end

struct writerContext
    aimWriterNum::Cint
    aimWriterName::NTuple{16, Ptr{Cchar}}
    aimWriterDLL::NTuple{16, Ptr{Cvoid}}
    aimExtension::NTuple{16, AIMext}
    aimWriter::NTuple{16, AIMwriter}
end

struct aimInfo
    magicnumber::Cint
    instance::Cint
    funID::Cint
    pIndex::Cint
    irow::Cint
    icol::Cint
    problem::Ptr{capsProblem}
    analysis::Ptr{Cvoid}
    errs::capsErrs
    wCntxt::writerContext
end

struct capsAnalysis
    loadName::Ptr{Cchar}
    fullPath::Ptr{Cchar}
    path::Ptr{Cchar}
    unitSys::Ptr{Cchar}
    major::Cint
    minor::Cint
    instStore::Ptr{Cvoid}
    autoexec::Cint
    eFlag::Cint
    reload::Cint
    intents::Ptr{Cchar}
    info::aimInfo
    nField::Cint
    fields::Ptr{Ptr{Cchar}}
    ranks::Ptr{Cint}
    fInOut::Ptr{Cint}
    nAnalysisIn::Cint
    analysisIn::Ptr{Ptr{capsObject}}
    nAnalysisOut::Cint
    analysisOut::Ptr{Ptr{capsObject}}
    nAnalysisDynO::Cint
    analysisDynO::Ptr{Ptr{capsObject}}
    nBody::Cint
    bodies::Ptr{ego}
    nTess::Cint
    tess::Ptr{ego}
    uSsN::Culong
    pre::capsOwn
end

struct capsBound
    dim::Cint
    state::Cint
    lunits::Ptr{Cchar}
    plimits::NTuple{4, Cdouble}
    geom::ego
    iBody::Cint
    iEnt::Cint
    curve::Ptr{capsAprx1D}
    surface::Ptr{capsAprx2D}
    index::Cint
    nVertexSet::Cint
    vertexSet::Ptr{Ptr{capsObject}}
end

struct capsVertexSet
    analysis::Ptr{capsObject}
    discr::Ptr{capsDiscr}
    nDataSets::Cint
    dataSets::Ptr{Ptr{capsObject}}
end

struct capsDataSet
    ftype::Cint
    npts::Cint
    rank::Cint
    data::Ptr{Cdouble}
    units::Ptr{Cchar}
    startup::Ptr{Cdouble}
    linkMethod::Cint
    link::Ptr{capsObject}
end

function EG_alloc(nbytes)
    ccall((:EG_alloc, "libegads"), Ptr{Cvoid}, (Cint,), nbytes)
end

function EG_calloc(nele, size)
    ccall((:EG_calloc, "libegads"), Ptr{Cvoid}, (Cint, Cint), nele, size)
end

function EG_reall(ptr, nbytes)
    ccall((:EG_reall, "libegads"), Ptr{Cvoid}, (Ptr{Cvoid}, Cint), ptr, nbytes)
end

function EG_strdup(str)
    ccall((:EG_strdup, "libegads"), Ptr{Cchar}, (Ptr{Cchar},), str)
end

function EG_free(pointer)
    ccall((:EG_free, "libegads"), Cvoid, (Ptr{Cvoid},), pointer)
end

function EG_revision(major, minor, OCCrev)
    ccall((:EG_revision, "libegads"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cchar}}), major, minor, OCCrev)
end

function EG_open(context)
    ccall((:EG_open, "libegads"), Cint, (Ptr{ego},), context)
end

function EG_loadModel(context, bflg, name, model)
    ccall((:EG_loadModel, "libegads"), Cint, (ego, Cint, Ptr{Cchar}, Ptr{ego}), context, bflg, name, model)
end

function EG_saveModel(model, name)
    ccall((:EG_saveModel, "libegads"), Cint, (ego, Ptr{Cchar}), model, name)
end

function EG_exportModel(model, nbytes, stream)
    ccall((:EG_exportModel, "libegads"), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cchar}}), model, nbytes, stream)
end

function EG_importModel(context, nbytes, stream, model)
    ccall((:EG_importModel, "libegads"), Cint, (ego, Cint, Ptr{Cchar}, Ptr{ego}), context, nbytes, stream, model)
end

function EG_deleteObject(object)
    ccall((:EG_deleteObject, "libegads"), Cint, (ego,), object)
end

function EG_makeTransform(context, xform, oform)
    ccall((:EG_makeTransform, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{ego}), context, xform, oform)
end

function EG_getTransformation(oform, xform)
    ccall((:EG_getTransformation, "libegads"), Cint, (ego, Ptr{Cdouble}), oform, xform)
end

function EG_getContext(object, context)
    ccall((:EG_getContext, "libegads"), Cint, (ego, Ptr{ego}), object, context)
end

function EG_setOutLevel(context, outLevel)
    ccall((:EG_setOutLevel, "libegads"), Cint, (ego, Cint), context, outLevel)
end

function EG_updateThread(context)
    ccall((:EG_updateThread, "libegads"), Cint, (ego,), context)
end

function EG_getInfo(object, oclass, mtype, topObj, prev, next)
    ccall((:EG_getInfo, "libegads"), Cint, (ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{ego}, Ptr{ego}), object, oclass, mtype, topObj, prev, next)
end

function EG_copyObject(object, oform, copy)
    ccall((:EG_copyObject, "libegads"), Cint, (ego, Ptr{Cvoid}, Ptr{ego}), object, oform, copy)
end

function EG_flipObject(object, flippedCopy)
    ccall((:EG_flipObject, "libegads"), Cint, (ego, Ptr{ego}), object, flippedCopy)
end

function EG_close(context)
    ccall((:EG_close, "libegads"), Cint, (ego,), context)
end

function EG_setUserPointer(context, ptr)
    ccall((:EG_setUserPointer, "libegads"), Cint, (ego, Ptr{Cvoid}), context, ptr)
end

function EG_getUserPointer(context, ptr)
    ccall((:EG_getUserPointer, "libegads"), Cint, (ego, Ptr{Ptr{Cvoid}}), context, ptr)
end

function EG_setFullAttrs(context, flag)
    ccall((:EG_setFullAttrs, "libegads"), Cint, (ego, Cint), context, flag)
end

function EG_attributeAdd(obj, name, type, len, ints, reals, str)
    ccall((:EG_attributeAdd, "libegads"), Cint, (ego, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cchar}), obj, name, type, len, ints, reals, str)
end

function EG_attributeDel(object, name)
    ccall((:EG_attributeDel, "libegads"), Cint, (ego, Ptr{Cchar}), object, name)
end

function EG_attributeNum(obj, num)
    ccall((:EG_attributeNum, "libegads"), Cint, (ego, Ptr{Cint}), obj, num)
end

function EG_attributeGet(obj, index, name, atype, len, ints, reals, str)
    ccall((:EG_attributeGet, "libegads"), Cint, (ego, Cint, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cchar}}), obj, index, name, atype, len, ints, reals, str)
end

function EG_attributeRet(obj, name, atype, len, ints, reals, str)
    ccall((:EG_attributeRet, "libegads"), Cint, (ego, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cchar}}), obj, name, atype, len, ints, reals, str)
end

function EG_attributeDup(src, dst)
    ccall((:EG_attributeDup, "libegads"), Cint, (ego, ego), src, dst)
end

function EG_attributeAddSeq(obj, name, type, len, ints, reals, str)
    ccall((:EG_attributeAddSeq, "libegads"), Cint, (ego, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cchar}), obj, name, type, len, ints, reals, str)
end

function EG_attributeNumSeq(obj, name, num)
    ccall((:EG_attributeNumSeq, "libegads"), Cint, (ego, Ptr{Cchar}, Ptr{Cint}), obj, name, num)
end

function EG_attributeRetSeq(obj, name, index, atype, len, ints, reals, str)
    ccall((:EG_attributeRetSeq, "libegads"), Cint, (ego, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cchar}}), obj, name, index, atype, len, ints, reals, str)
end

function EG_getGeometry(geom, oclass, mtype, refGeom, ivec, rvec)
    ccall((:EG_getGeometry, "libegads"), Cint, (ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), geom, oclass, mtype, refGeom, ivec, rvec)
end

function EG_makeGeometry(context, oclass, mtype, refGeom, ivec, rvec, geom)
    ccall((:EG_makeGeometry, "libegads"), Cint, (ego, Cint, Cint, ego, Ptr{Cint}, Ptr{Cdouble}, Ptr{ego}), context, oclass, mtype, refGeom, ivec, rvec, geom)
end

function EG_getRange(geom, range, periodic)
    ccall((:EG_getRange, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{Cint}), geom, range, periodic)
end

function EG_evaluate(geom, param, results)
    ccall((:EG_evaluate, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}), geom, param, results)
end

function EG_invEvaluate(geom, xyz, param, results)
    ccall((:EG_invEvaluate, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), geom, xyz, param, results)
end

function EG_invEvaluateGuess(geom, xyz, param, results)
    ccall((:EG_invEvaluateGuess, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), geom, xyz, param, results)
end

function EG_arcLength(geom, t1, t2, alen)
    ccall((:EG_arcLength, "libegads"), Cint, (ego, Cdouble, Cdouble, Ptr{Cdouble}), geom, t1, t2, alen)
end

function EG_curvature(geom, param, results)
    ccall((:EG_curvature, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}), geom, param, results)
end

function EG_approximate(context, maxdeg, tol, sizes, xyzs, bspline)
    ccall((:EG_approximate, "libegads"), Cint, (ego, Cint, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Ptr{ego}), context, maxdeg, tol, sizes, xyzs, bspline)
end

function EG_fitTriangles(context, npts, xyzs, ntris, tris, tric, tol, bspline)
    ccall((:EG_fitTriangles, "libegads"), Cint, (ego, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Cdouble, Ptr{ego}), context, npts, xyzs, ntris, tris, tric, tol, bspline)
end

function EG_otherCurve(surface, curve, tol, newcurve)
    ccall((:EG_otherCurve, "libegads"), Cint, (ego, ego, Cdouble, Ptr{ego}), surface, curve, tol, newcurve)
end

function EG_isSame(geom1, geom2)
    ccall((:EG_isSame, "libegads"), Cint, (ego, ego), geom1, geom2)
end

function EG_isoCline(surface, UV, val, newcurve)
    ccall((:EG_isoCline, "libegads"), Cint, (ego, Cint, Cdouble, Ptr{ego}), surface, UV, val, newcurve)
end

function EG_convertToBSpline(geom, bspline)
    ccall((:EG_convertToBSpline, "libegads"), Cint, (ego, Ptr{ego}), geom, bspline)
end

function EG_convertToBSplineRange(geom, range, bspline)
    ccall((:EG_convertToBSplineRange, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{ego}), geom, range, bspline)
end

function EG_skinning(nCurves, curves, skinning_degree, surface)
    ccall((:EG_skinning, "libegads"), Cint, (Cint, Ptr{ego}, Cint, Ptr{ego}), nCurves, curves, skinning_degree, surface)
end

function EG_tolerance(topo, tol)
    ccall((:EG_tolerance, "libegads"), Cint, (ego, Ptr{Cdouble}), topo, tol)
end

function EG_getTolerance(topo, tol)
    ccall((:EG_getTolerance, "libegads"), Cint, (ego, Ptr{Cdouble}), topo, tol)
end

function EG_setTolerance(topo, tol)
    ccall((:EG_setTolerance, "libegads"), Cint, (ego, Cdouble), topo, tol)
end

function EG_getTopology(topo, geom, oclass, type, limits, nChildren, children, sense)
    ccall((:EG_getTopology, "libegads"), Cint, (ego, Ptr{ego}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Ptr{Cint}}), topo, geom, oclass, type, limits, nChildren, children, sense)
end

function EG_makeTopology(context, geom, oclass, mtype, limits, nChildren, children, senses, topo)
    ccall((:EG_makeTopology, "libegads"), Cint, (ego, ego, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{ego}, Ptr{Cint}, Ptr{ego}), context, geom, oclass, mtype, limits, nChildren, children, senses, topo)
end

function EG_makeLoop(nedge, edges, geom, toler, result)
    ccall((:EG_makeLoop, "libegads"), Cint, (Cint, Ptr{ego}, ego, Cdouble, Ptr{ego}), nedge, edges, geom, toler, result)
end

function EG_getArea(object, limits, area)
    ccall((:EG_getArea, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}), object, limits, area)
end

function EG_makeFace(object, mtype, limits, face)
    ccall((:EG_makeFace, "libegads"), Cint, (ego, Cint, Ptr{Cdouble}, Ptr{ego}), object, mtype, limits, face)
end

function EG_getBodyTopos(body, src, oclass, ntopo, topos)
    ccall((:EG_getBodyTopos, "libegads"), Cint, (ego, ego, Cint, Ptr{Cint}, Ptr{Ptr{ego}}), body, src, oclass, ntopo, topos)
end

function EG_indexBodyTopo(body, src)
    ccall((:EG_indexBodyTopo, "libegads"), Cint, (ego, ego), body, src)
end

function EG_objectBodyTopo(body, oclass, index, obj)
    ccall((:EG_objectBodyTopo, "libegads"), Cint, (ego, Cint, Cint, Ptr{ego}), body, oclass, index, obj)
end

function EG_inTopology(topo, xyz)
    ccall((:EG_inTopology, "libegads"), Cint, (ego, Ptr{Cdouble}), topo, xyz)
end

function EG_inFace(face, uv)
    ccall((:EG_inFace, "libegads"), Cint, (ego, Ptr{Cdouble}), face, uv)
end

function EG_getEdgeUV(face, edge, sense, t, UV)
    ccall((:EG_getEdgeUV, "libegads"), Cint, (ego, ego, Cint, Cdouble, Ptr{Cdouble}), face, edge, sense, t, UV)
end

function EG_getEdgeUVs(face, edge, sense, nts, ts, UV)
    ccall((:EG_getEdgeUVs, "libegads"), Cint, (ego, ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), face, edge, sense, nts, ts, UV)
end

function EG_getPCurve(face, edge, sense, mtype, ivec, rvec)
    ccall((:EG_getPCurve, "libegads"), Cint, (ego, ego, Cint, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), face, edge, sense, mtype, ivec, rvec)
end

function EG_getWindingAngle(edge, t, angle)
    ccall((:EG_getWindingAngle, "libegads"), Cint, (ego, Cdouble, Ptr{Cdouble}), edge, t, angle)
end

function EG_getBody(object, body)
    ccall((:EG_getBody, "libegads"), Cint, (ego, Ptr{ego}), object, body)
end

function EG_makeSolidBody(context, stype, rvec, body)
    ccall((:EG_makeSolidBody, "libegads"), Cint, (ego, Cint, Ptr{Cdouble}, Ptr{ego}), context, stype, rvec, body)
end

function EG_getBoundingBox(topo, bbox)
    ccall((:EG_getBoundingBox, "libegads"), Cint, (ego, Ptr{Cdouble}), topo, bbox)
end

function EG_getMassProperties(topo, result)
    ccall((:EG_getMassProperties, "libegads"), Cint, (ego, Ptr{Cdouble}), topo, result)
end

function EG_isEquivalent(topo1, topo2)
    ccall((:EG_isEquivalent, "libegads"), Cint, (ego, ego), topo1, topo2)
end

function EG_sewFaces(nobj, objs, toler, flag, result)
    ccall((:EG_sewFaces, "libegads"), Cint, (Cint, Ptr{ego}, Cdouble, Cint, Ptr{ego}), nobj, objs, toler, flag, result)
end

function EG_makeNmWireBody(nobj, objs, toler, result)
    ccall((:EG_makeNmWireBody, "libegads"), Cint, (Cint, Ptr{ego}, Cdouble, Ptr{ego}), nobj, objs, toler, result)
end

function EG_replaceFaces(body, nobj, objs, result)
    ccall((:EG_replaceFaces, "libegads"), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}), body, nobj, objs, result)
end

function EG_mapBody(sBody, dBody, fAttr, mapBody)
    ccall((:EG_mapBody, "libegads"), Cint, (ego, ego, Ptr{Cchar}, Ptr{ego}), sBody, dBody, fAttr, mapBody)
end

function EG_matchBodyEdges(body1, body2, toler, nmatch, match)
    ccall((:EG_matchBodyEdges, "libegads"), Cint, (ego, ego, Cdouble, Ptr{Cint}, Ptr{Ptr{Cint}}), body1, body2, toler, nmatch, match)
end

function EG_matchBodyFaces(body1, body2, toler, nmatch, match)
    ccall((:EG_matchBodyFaces, "libegads"), Cint, (ego, ego, Cdouble, Ptr{Cint}, Ptr{Ptr{Cint}}), body1, body2, toler, nmatch, match)
end

function EG_setTessParam(context, iparam, value, oldvalue)
    ccall((:EG_setTessParam, "libegads"), Cint, (ego, Cint, Cdouble, Ptr{Cdouble}), context, iparam, value, oldvalue)
end

function EG_makeTessGeom(obj, params, sizes, tess)
    ccall((:EG_makeTessGeom, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{Cint}, Ptr{ego}), obj, params, sizes, tess)
end

function EG_getTessGeom(tess, sizes, xyz)
    ccall((:EG_getTessGeom, "libegads"), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cdouble}}), tess, sizes, xyz)
end

function EG_makeTessBody(object, params, tess)
    ccall((:EG_makeTessBody, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{ego}), object, params, tess)
end

function EG_remakeTess(tess, nobj, objs, params)
    ccall((:EG_remakeTess, "libegads"), Cint, (ego, Cint, Ptr{ego}, Ptr{Cdouble}), tess, nobj, objs, params)
end

function EG_finishTess(tess, params)
    ccall((:EG_finishTess, "libegads"), Cint, (ego, Ptr{Cdouble}), tess, params)
end

function EG_mapTessBody(tess, body, mapTess)
    ccall((:EG_mapTessBody, "libegads"), Cint, (ego, ego, Ptr{ego}), tess, body, mapTess)
end

function EG_locateTessBody(tess, npt, ifaces, uv, itri, results)
    ccall((:EG_locateTessBody, "libegads"), Cint, (ego, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}), tess, npt, ifaces, uv, itri, results)
end

function EG_getTessEdge(tess, eIndex, len, xyz, t)
    ccall((:EG_getTessEdge, "libegads"), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}), tess, eIndex, len, xyz, t)
end

function EG_getTessFace(tess, fIndex, len, xyz, uv, ptype, pindex, ntri, tris, tric)
    ccall((:EG_getTessFace, "libegads"), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}), tess, fIndex, len, xyz, uv, ptype, pindex, ntri, tris, tric)
end

function EG_getTessLoops(tess, fIndex, nloop, lIndices)
    ccall((:EG_getTessLoops, "libegads"), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cint}}), tess, fIndex, nloop, lIndices)
end

function EG_getTessQuads(tess, nquad, fIndices)
    ccall((:EG_getTessQuads, "libegads"), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cint}}), tess, nquad, fIndices)
end

function EG_makeQuads(tess, params, fIndex)
    ccall((:EG_makeQuads, "libegads"), Cint, (ego, Ptr{Cdouble}, Cint), tess, params, fIndex)
end

function EG_getQuads(tess, fIndex, len, xyz, uv, ptype, pindex, npatch)
    ccall((:EG_getQuads, "libegads"), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}), tess, fIndex, len, xyz, uv, ptype, pindex, npatch)
end

function EG_getPatch(tess, fIndex, patch, nu, nv, ipts, bounds)
    ccall((:EG_getPatch, "libegads"), Cint, (ego, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}), tess, fIndex, patch, nu, nv, ipts, bounds)
end

function EG_quadTess(tess, quadTess)
    ccall((:EG_quadTess, "libegads"), Cint, (ego, Ptr{ego}), tess, quadTess)
end

function EG_insertEdgeVerts(tess, eIndex, vIndex, npts, t)
    ccall((:EG_insertEdgeVerts, "libegads"), Cint, (ego, Cint, Cint, Cint, Ptr{Cdouble}), tess, eIndex, vIndex, npts, t)
end

function EG_deleteEdgeVert(tess, eIndex, vIndex, dir)
    ccall((:EG_deleteEdgeVert, "libegads"), Cint, (ego, Cint, Cint, Cint), tess, eIndex, vIndex, dir)
end

function EG_moveEdgeVert(tess, eIndex, vIndex, t)
    ccall((:EG_moveEdgeVert, "libegads"), Cint, (ego, Cint, Cint, Cdouble), tess, eIndex, vIndex, t)
end

function EG_openTessBody(tess)
    ccall((:EG_openTessBody, "libegads"), Cint, (ego,), tess)
end

function EG_initTessBody(object, tess)
    ccall((:EG_initTessBody, "libegads"), Cint, (ego, Ptr{ego}), object, tess)
end

function EG_statusTessBody(tess, body, state, np)
    ccall((:EG_statusTessBody, "libegads"), Cint, (ego, Ptr{ego}, Ptr{Cint}, Ptr{Cint}), tess, body, state, np)
end

function EG_setTessEdge(tess, eIndex, len, xyz, t)
    ccall((:EG_setTessEdge, "libegads"), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), tess, eIndex, len, xyz, t)
end

function EG_setTessFace(tess, fIndex, len, xyz, uv, ntri, tris)
    ccall((:EG_setTessFace, "libegads"), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}), tess, fIndex, len, xyz, uv, ntri, tris)
end

function EG_localToGlobal(tess, index, _local, _global)
    ccall((:EG_localToGlobal, "libegads"), Cint, (ego, Cint, Cint, Ptr{Cint}), tess, index, _local, _global)
end

function EG_getGlobal(tess, _global, ptype, pindex, xyz)
    ccall((:EG_getGlobal, "libegads"), Cint, (ego, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), tess, _global, ptype, pindex, xyz)
end

function EG_saveTess(tess, name)
    ccall((:EG_saveTess, "libegads"), Cint, (ego, Ptr{Cchar}), tess, name)
end

function EG_loadTess(body, name, tess)
    ccall((:EG_loadTess, "libegads"), Cint, (ego, Ptr{Cchar}, Ptr{ego}), body, name, tess)
end

function EG_tessMassProps(tess, props)
    ccall((:EG_tessMassProps, "libegads"), Cint, (ego, Ptr{Cdouble}), tess, props)
end

function EG_fuseSheets(src, tool, sheet)
    ccall((:EG_fuseSheets, "libegads"), Cint, (ego, ego, Ptr{ego}), src, tool, sheet)
end

function EG_generalBoolean(src, tool, oper, tol, model)
    ccall((:EG_generalBoolean, "libegads"), Cint, (ego, ego, Cint, Cdouble, Ptr{ego}), src, tool, oper, tol, model)
end

function EG_solidBoolean(src, tool, oper, model)
    ccall((:EG_solidBoolean, "libegads"), Cint, (ego, ego, Cint, Ptr{ego}), src, tool, oper, model)
end

function EG_intersection(src, tool, nedge, facEdg, model)
    ccall((:EG_intersection, "libegads"), Cint, (ego, ego, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{ego}), src, tool, nedge, facEdg, model)
end

function EG_imprintBody(src, nedge, facEdg, result)
    ccall((:EG_imprintBody, "libegads"), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}), src, nedge, facEdg, result)
end

function EG_filletBody(src, nedge, edges, radius, result, facemap)
    ccall((:EG_filletBody, "libegads"), Cint, (ego, Cint, Ptr{ego}, Cdouble, Ptr{ego}, Ptr{Ptr{Cint}}), src, nedge, edges, radius, result, facemap)
end

function EG_chamferBody(src, nedge, edges, faces, dis1, dis2, result, facemap)
    ccall((:EG_chamferBody, "libegads"), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}, Cdouble, Cdouble, Ptr{ego}, Ptr{Ptr{Cint}}), src, nedge, edges, faces, dis1, dis2, result, facemap)
end

function EG_hollowBody(src, nface, faces, offset, join, result, facemap)
    ccall((:EG_hollowBody, "libegads"), Cint, (ego, Cint, Ptr{ego}, Cdouble, Cint, Ptr{ego}, Ptr{Ptr{Cint}}), src, nface, faces, offset, join, result, facemap)
end

function EG_extrude(src, dist, dir, result)
    ccall((:EG_extrude, "libegads"), Cint, (ego, Cdouble, Ptr{Cdouble}, Ptr{ego}), src, dist, dir, result)
end

function EG_rotate(src, angle, axis, result)
    ccall((:EG_rotate, "libegads"), Cint, (ego, Cdouble, Ptr{Cdouble}, Ptr{ego}), src, angle, axis, result)
end

function EG_sweep(src, spine, mode, result)
    ccall((:EG_sweep, "libegads"), Cint, (ego, ego, Cint, Ptr{ego}), src, spine, mode, result)
end

function EG_loft(nsec, secs, opt, result)
    ccall((:EG_loft, "libegads"), Cint, (Cint, Ptr{ego}, Cint, Ptr{ego}), nsec, secs, opt, result)
end

function EG_blend(nsec, secs, rc1, rcN, result)
    ccall((:EG_blend, "libegads"), Cint, (Cint, Ptr{ego}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ego}), nsec, secs, rc1, rcN, result)
end

function EG_ruled(nsec, secs, result)
    ccall((:EG_ruled, "libegads"), Cint, (Cint, Ptr{ego}, Ptr{ego}), nsec, secs, result)
end

function EG_initEBody(tess, angle, ebody)
    ccall((:EG_initEBody, "libegads"), Cint, (ego, Cdouble, Ptr{ego}), tess, angle, ebody)
end

function EG_finishEBody(EBody)
    ccall((:EG_finishEBody, "libegads"), Cint, (ego,), EBody)
end

function EG_makeEFace(EBody, nFace, Faces, EFace)
    ccall((:EG_makeEFace, "libegads"), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}), EBody, nFace, Faces, EFace)
end

function EG_makeAttrEFaces(EBody, attrName, nEFace, EFaces)
    ccall((:EG_makeAttrEFaces, "libegads"), Cint, (ego, Ptr{Cchar}, Ptr{Cint}, Ptr{Ptr{ego}}), EBody, attrName, nEFace, EFaces)
end

function EG_effectiveMap(EObject, eparam, Object, param)
    ccall((:EG_effectiveMap, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{ego}, Ptr{Cdouble}), EObject, eparam, Object, param)
end

function EG_effectiveEdgeList(EEdge, nedge, edges, senses, tstart)
    ccall((:EG_effectiveEdgeList, "libegads"), Cint, (ego, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), EEdge, nedge, edges, senses, tstart)
end

function EG_effectiveTri(EObj, uv, fIndex, itri, w)
    ccall((:EG_effectiveTri, "libegads"), Cint, (ego, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), EObj, uv, fIndex, itri, w)
end

function EG_inTriExact(t1, t2, t3, p, w)
    ccall((:EG_inTriExact, "libegads"), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), t1, t2, t3, p, w)
end

struct attr_T
    name::Ptr{Cchar}
    defn::Ptr{Cchar}
    type::Cint
end

struct grat_T
    object::Ptr{Cvoid}
    active::Cint
    color::Cint
    bcolor::Cint
    mcolor::Cint
    lwidth::Cint
    ptsize::Cint
    render::Cint
    dirty::Cint
end

struct varg_T
    data::NTuple{32, UInt8}
end

function Base.getproperty(x::Ptr{varg_T}, f::Symbol)
    f === :nval && return Ptr{Cint}(x + 0)
    f === :nrow && return Ptr{Cint}(x + 4)
    f === :ncol && return Ptr{Cint}(x + 8)
    f === :str && return Ptr{Ptr{Cchar}}(x + 16)
    f === :val && return Ptr{Ptr{Cdouble}}(x + 16)
    f === :dot && return Ptr{Ptr{Cdouble}}(x + 24)
    return getfield(x, f)
end

function Base.getproperty(x::varg_T, f::Symbol)
    r = Ref{varg_T}(x)
    ptr = Base.unsafe_convert(Ptr{varg_T}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{varg_T}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct mprp_T
    name::NTuple{10, Cchar}
    val::Cdouble
end

struct node_T
    nedge::Cint
    x::Cdouble
    y::Cdouble
    z::Cdouble
    ibody::Cint
    gratt::grat_T
    dxyz::Ptr{Cdouble}
    enode::ego
end

struct edge_T
    itype::Cint
    ibeg::Cint
    iend::Cint
    ileft::Cint
    irite::Cint
    nface::Cint
    ibody::Cint
    iford::Cint
    imark::Cint
    gratt::grat_T
    dxyz::Ptr{Cdouble}
    dt::Ptr{Cdouble}
    globid::Cint
    eedge::ego
end

struct face_T
    ibody::Cint
    iford::Cint
    imark::Cint
    gratt::grat_T
    eggdata::Ptr{Cvoid}
    dxyz::Ptr{Cdouble}
    duv::Ptr{Cdouble}
    globid::Cint
    eface::ego
end

struct body_T
    ibrch::Cint
    brtype::Cint
    ileft::Cint
    irite::Cint
    ichld::Cint
    igroup::Cint
    arg::NTuple{10, varg_T}
    ebody::ego
    eebody::ego
    etess::ego
    eetess::ego
    npnts::Cint
    ntris::Cint
    onstack::Cint
    hasdots::Cint
    hasdxyz::Cint
    botype::Cint
    nonmani::Cint
    CPU::Cdouble
    nnode::Cint
    node::Ptr{node_T}
    nedge::Cint
    edge::Ptr{edge_T}
    nface::Cint
    face::Ptr{face_T}
    hassens::Cint
    gratt::grat_T
end

struct brch_T
    name::Ptr{Cchar}
    type::Cint
    bclass::Cint
    level::Cint
    indent::Cint
    filename::Ptr{Cchar}
    linenum::Cint
    actv::Cint
    dirty::Cint
    nattr::Cint
    attr::Ptr{attr_T}
    ileft::Cint
    irite::Cint
    ichld::Cint
    nmprp::Cint
    mprp::Ptr{mprp_T}
    narg::Cint
    arg1::Ptr{Cchar}
    arg2::Ptr{Cchar}
    arg3::Ptr{Cchar}
    arg4::Ptr{Cchar}
    arg5::Ptr{Cchar}
    arg6::Ptr{Cchar}
    arg7::Ptr{Cchar}
    arg8::Ptr{Cchar}
    arg9::Ptr{Cchar}
end

struct pmtr_T
    name::Ptr{Cchar}
    type::Cint
    scope::Cint
    mprop::Cint
    nrow::Cint
    ncol::Cint
    value::Ptr{Cdouble}
    dot::Ptr{Cdouble}
    lbnd::Ptr{Cdouble}
    ubnd::Ptr{Cdouble}
    str::Ptr{Cchar}
end

struct stor_T
    name::NTuple{65, Cchar}
    index::Cint
    nbody::Cint
    ibody::Ptr{Cint}
    ebody::Ptr{ego}
end

struct prof_T
    ncall::Cint
    time::Cint
end

# typedef void ( * mesgCB_H ) ( char message [ ] )
const mesgCB_H = Ptr{Cvoid}

# typedef void ( * sizeCB_H ) ( void * modl , int ipmtr , int nrow , int ncol )
const sizeCB_H = Ptr{Cvoid}

# typedef int ( * eggGenerate_H ) ( double [ ] , int [ ] , void * * )
const eggGenerate_H = Ptr{Cvoid}

# typedef int ( * eggMorph_H ) ( void * , double * , void * * )
const eggMorph_H = Ptr{Cvoid}

# typedef int ( * eggInfo_H ) ( void * , int * , int * , const double * * , const int * * , int * , const int * * )
const eggInfo_H = Ptr{Cvoid}

# typedef int ( * eggDump_H ) ( void * , FILE * )
const eggDump_H = Ptr{Cvoid}

# typedef int ( * eggLoad_H ) ( FILE * , void * * )
const eggLoad_H = Ptr{Cvoid}

# typedef int ( * eggFree_H ) ( void * )
const eggFree_H = Ptr{Cvoid}

mutable struct modl_T
    magic::Cint
    checked::Cint
    embedded::Cint
    ibrch::Cint
    nextseq::Cint
    ngroup::Cint
    recycle::Cint
    verify::Cint
    cleanup::Cint
    dumpEgads::Cint
    loadEgads::Cint
    hasMPs::Cint
    printStack::Cint
    tessAtEnd::Cint
    erepAtEnd::Cint
    bodyLoaded::Cint
    seltype::Cint
    selbody::Cint
    selsize::Cint
    sellist::Ptr{Cint}
    level::Cint
    scope::NTuple{11, Cint}
    filename::NTuple{256, Cchar}
    tmpDirNum::Cint
    ninline::Cint
    minline::Cint
    sinline::Ptr{Cchar}
    nattr::Cint
    attr::Ptr{attr_T}
    nstor::Cint
    stor::Ptr{stor_T}
    nbrch::Cint
    mbrch::Cint
    brch::Ptr{brch_T}
    npmtr::Cint
    mpmtr::Cint
    pmtr::Ptr{pmtr_T}
    nbody::Cint
    mbody::Cint
    body::Ptr{body_T}
    forceFDs::Cdouble
    needFDs::Cint
    numdots::Cint
    perturb::Ptr{modl_T}
    basemodl::Ptr{modl_T}
    dtime::Cdouble
    context::ego
    userdata::Ptr{Cvoid}
    mesgCB::mesgCB_H
    bcstCB::mesgCB_H
    sizeCB::sizeCB_H
    eggname::NTuple{256, Cchar}
    eggGenerate::eggGenerate_H
    eggMorph::eggMorph_H
    eggInfo::eggInfo_H
    eggDump::eggDump_H
    eggLoad::eggLoad_H
    eggFree::eggFree_H
    nwarn::Cint
    sigCode::Cint
    sigMesg::Ptr{Cchar}
    profile::NTuple{101, prof_T}
end

function ocsmVersion(imajor, iminor)
    ccall((:ocsmVersion, "libocsm"), Cint, (Ptr{Cint}, Ptr{Cint}), imajor, iminor)
end

function ocsmSetOutLevel(ilevel)
    ccall((:ocsmSetOutLevel, "libocsm"), Cint, (Cint,), ilevel)
end

function ocsmSetAuxPtr(newAuxPtr)
    ccall((:ocsmSetAuxPtr, "libocsm"), Cint, (Ptr{Cvoid},), newAuxPtr)
end

function ocsmGetAuxPtr(oldAuxPtr)
    ccall((:ocsmGetAuxPtr, "libocsm"), Cint, (Ptr{Ptr{Cvoid}},), oldAuxPtr)
end

function ocsmLoad(filename, modl)
    ccall((:ocsmLoad, "libocsm"), Cint, (Ptr{Cchar}, Ptr{Ptr{Cvoid}}), filename, modl)
end

function ocsmLoadFromModel(emodel, modl)
    ccall((:ocsmLoadFromModel, "libocsm"), Cint, (ego, Ptr{Ptr{Cvoid}}), emodel, modl)
end

function ocsmLoadDict(modl, dictname)
    ccall((:ocsmLoadDict, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, dictname)
end

function ocsmAdjustUDCs(modl)
    ccall((:ocsmAdjustUDCs, "libocsm"), Cint, (Ptr{Cvoid},), modl)
end

function ocsmUpdateDespmtrs(modl, filename)
    ccall((:ocsmUpdateDespmtrs, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, filename)
end

function ocsmGetFilelist(modl, filelist)
    ccall((:ocsmGetFilelist, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Ptr{Cchar}}), modl, filelist)
end

function ocsmGetFiletree(modl, info)
    ccall((:ocsmGetFiletree, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Ptr{Cchar}}), modl, info)
end

function ocsmSave(modl, filename)
    ccall((:ocsmSave, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, filename)
end

function ocsmSaveDespmtrs(modl, filename)
    ccall((:ocsmSaveDespmtrs, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, filename)
end

function ocsmCopy(srcModl, newModl)
    ccall((:ocsmCopy, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Ptr{Cvoid}}), srcModl, newModl)
end

function ocsmFree(modl)
    ccall((:ocsmFree, "libocsm"), Cint, (Ptr{Cvoid},), modl)
end

function ocsmInfo(modl, nbrch, npmtr, nbody)
    ccall((:ocsmInfo, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), modl, nbrch, npmtr, nbody)
end

function ocsmCheck(modl)
    ccall((:ocsmCheck, "libocsm"), Cint, (Ptr{Cvoid},), modl)
end

function ocsmRegMesgCB(modl, callback)
    ccall((:ocsmRegMesgCB, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), modl, callback)
end

function ocsmRegBcstCB(modl, callback)
    ccall((:ocsmRegBcstCB, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), modl, callback)
end

function ocsmRegSizeCB(modl, callback)
    ccall((:ocsmRegSizeCB, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), modl, callback)
end

function ocsmBuild(modl, buildTo, builtTo, nbody, body)
    ccall((:ocsmBuild, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), modl, buildTo, builtTo, nbody, body)
end

function ocsmPrintProfile(modl, filename)
    ccall((:ocsmPrintProfile, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, filename)
end

function ocsmTessellate(modl, ibody)
    ccall((:ocsmTessellate, "libocsm"), Cint, (Ptr{Cvoid}, Cint), modl, ibody)
end

function ocsmBodyDetails(modl, fiename, linenum, info)
    ccall((:ocsmBodyDetails, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Ptr{Cchar}}), modl, fiename, linenum, info)
end

function ocsmPerturb(modl, npmtrs, ipmtrs, irows, icols, values)
    ccall((:ocsmPerturb, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), modl, npmtrs, ipmtrs, irows, icols, values)
end

function ocsmUpdateTess(modl, ibody, filename)
    ccall((:ocsmUpdateTess, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}), modl, ibody, filename)
end

function ocsmClearance(modl, ibody1, ibody2, dist, pnt1, pnt2)
    ccall((:ocsmClearance, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), modl, ibody1, ibody2, dist, pnt1, pnt2)
end

function ocsmShowTblOfContents(modl, irow, icol, info)
    ccall((:ocsmShowTblOfContents, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Ptr{Cchar}}), modl, irow, icol, info)
end

function ocsmNewBrch(modl, iafter, type, filename, linenum, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    ccall((:ocsmNewBrch, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cchar}, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}), modl, iafter, type, filename, linenum, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
end

function ocsmGetBrch(modl, ibrch, type, bclass, actv, ichld, ileft, irite, narg, nattr)
    ccall((:ocsmGetBrch, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), modl, ibrch, type, bclass, actv, ichld, ileft, irite, narg, nattr)
end

function ocsmSetBrch(modl, ibrch, actv)
    ccall((:ocsmSetBrch, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint), modl, ibrch, actv)
end

function ocsmDelBrch(modl, ibrch)
    ccall((:ocsmDelBrch, "libocsm"), Cint, (Ptr{Cvoid}, Cint), modl, ibrch)
end

function ocsmPrintBrchs(modl, filename)
    ccall((:ocsmPrintBrchs, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, filename)
end

function ocsmGetArg(modl, ibrch, iarg, defn, value, dot)
    ccall((:ocsmGetArg, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}), modl, ibrch, iarg, defn, value, dot)
end

function ocsmSetArg(modl, ibrch, iarg, defn)
    ccall((:ocsmSetArg, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cchar}), modl, ibrch, iarg, defn)
end

function ocsmRetAttr(modl, ibrch, iattr, aname, avalue)
    ccall((:ocsmRetAttr, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}), modl, ibrch, iattr, aname, avalue)
end

function ocsmGetAttr(modl, ibrch, aname, avalue)
    ccall((:ocsmGetAttr, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cchar}), modl, ibrch, aname, avalue)
end

function ocsmSetAttr(modl, ibrch, aname, avalue)
    ccall((:ocsmSetAttr, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cchar}), modl, ibrch, aname, avalue)
end

function ocsmTraceAttrs(modl, pattern, info)
    ccall((:ocsmTraceAttrs, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), modl, pattern, info)
end

function ocsmRetCsys(modl, ibrch, icsys, cname, cvalue)
    ccall((:ocsmRetCsys, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}), modl, ibrch, icsys, cname, cvalue)
end

function ocsmGetCsys(modl, ibrch, cname, cvalue)
    ccall((:ocsmGetCsys, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cchar}), modl, ibrch, cname, cvalue)
end

function ocsmSetCsys(modl, ibrch, cname, cvalue)
    ccall((:ocsmSetCsys, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cchar}), modl, ibrch, cname, cvalue)
end

function ocsmPrintAttrs(modl, filename)
    ccall((:ocsmPrintAttrs, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, filename)
end

function ocsmGetName(modl, ibrch, name)
    ccall((:ocsmGetName, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}), modl, ibrch, name)
end

function ocsmSetName(modl, ibrch, name)
    ccall((:ocsmSetName, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}), modl, ibrch, name)
end

function ocsmGetSketch(modl, ibrch, maxlen, begs, vars, cons, segs)
    ccall((:ocsmGetSketch, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}), modl, ibrch, maxlen, begs, vars, cons, segs)
end

function ocsmSolveSketch(modl, vars_in, cons, vars_out)
    ccall((:ocsmSolveSketch, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}), modl, vars_in, cons, vars_out)
end

function ocsmSaveSketch(modl, ibrch, vars, cons, segs)
    ccall((:ocsmSaveSketch, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}), modl, ibrch, vars, cons, segs)
end

function ocsmMakeEBody(modl, ibody, dihedral, entList)
    ccall((:ocsmMakeEBody, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cdouble, Ptr{Cchar}), modl, ibody, dihedral, entList)
end

function ocsmNewPmtr(modl, name, type, nrow, ncol)
    ccall((:ocsmNewPmtr, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Cint), modl, name, type, nrow, ncol)
end

function ocsmDelPmtr(modl, ipmtr)
    ccall((:ocsmDelPmtr, "libocsm"), Cint, (Ptr{Cvoid}, Cint), modl, ipmtr)
end

function ocsmFindPmtr(modl, name, type, nrow, ncol, ipmtr)
    ccall((:ocsmFindPmtr, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Cint, Ptr{Cint}), modl, name, type, nrow, ncol, ipmtr)
end

function ocsmGetPmtr(modl, ipmtr, type, nrow, ncol, name)
    ccall((:ocsmGetPmtr, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cchar}), modl, ipmtr, type, nrow, ncol, name)
end

function ocsmPrintPmtrs(modl, filename)
    ccall((:ocsmPrintPmtrs, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, filename)
end

function ocsmTracePmtrs(modl, pattern, info)
    ccall((:ocsmTracePmtrs, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), modl, pattern, info)
end

function ocsmGetValu(modl, ipmtr, irow, icol, value, dot)
    ccall((:ocsmGetValu, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), modl, ipmtr, irow, icol, value, dot)
end

function ocsmGetValuS(modl, ipmtr, str)
    ccall((:ocsmGetValuS, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}), modl, ipmtr, str)
end

function ocsmSetValu(modl, ipmtr, irow, icol, defn)
    ccall((:ocsmSetValu, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}), modl, ipmtr, irow, icol, defn)
end

function ocsmSetValuD(modl, ipmtr, irow, icol, value)
    ccall((:ocsmSetValuD, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Cdouble), modl, ipmtr, irow, icol, value)
end

function ocsmGetBnds(modl, ipmtr, irow, icol, lbound, ubound)
    ccall((:ocsmGetBnds, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), modl, ipmtr, irow, icol, lbound, ubound)
end

function ocsmSetBnds(modl, ipmtr, irow, icol, lbound, ubound)
    ccall((:ocsmSetBnds, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Cdouble, Cdouble), modl, ipmtr, irow, icol, lbound, ubound)
end

function ocsmSetDtime(modl, dtime)
    ccall((:ocsmSetDtime, "libocsm"), Cint, (Ptr{Cvoid}, Cdouble), modl, dtime)
end

function ocsmSetVel(modl, ipmtr, irow, icol, defn)
    ccall((:ocsmSetVel, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}), modl, ipmtr, irow, icol, defn)
end

function ocsmSetVelD(modl, ipmtr, irow, icol, dot)
    ccall((:ocsmSetVelD, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Cdouble), modl, ipmtr, irow, icol, dot)
end

function ocsmGetEgo(modl, ibody, seltype, iselect, theEgo)
    ccall((:ocsmGetEgo, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{ego}), modl, ibody, seltype, iselect, theEgo)
end

function ocsmSetEgo(modl, ibody, iselect, theEgo)
    ccall((:ocsmSetEgo, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, ego), modl, ibody, iselect, theEgo)
end

function ocsmFindEnt(modl, ibody, seltype, entID, ient)
    ccall((:ocsmFindEnt, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cint}, Ptr{Cint}), modl, ibody, seltype, entID, ient)
end

function ocsmGetUV(modl, ibody, seltype, iselect, npnt, xyz, uv)
    ccall((:ocsmGetUV, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), modl, ibody, seltype, iselect, npnt, xyz, uv)
end

function ocsmGetXYZ(modl, ibody, seltype, iselect, npnt, uv, xyz)
    ccall((:ocsmGetXYZ, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), modl, ibody, seltype, iselect, npnt, uv, xyz)
end

function ocsmGetNorm(modl, ibody, iface, npnt, uv, norm)
    ccall((:ocsmGetNorm, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), modl, ibody, iface, npnt, uv, norm)
end

function ocsmGetVel(modl, ibody, seltype, iselect, npnt, uv, vel)
    ccall((:ocsmGetVel, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), modl, ibody, seltype, iselect, npnt, uv, vel)
end

function ocsmSetEgg(modl, eggname)
    ccall((:ocsmSetEgg, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, eggname)
end

function ocsmGetTessNpnt(modl, ibody, seltype, iselect, npnt)
    ccall((:ocsmGetTessNpnt, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cint}), modl, ibody, seltype, iselect, npnt)
end

function ocsmGetTessVel(modl, ibody, seltype, iselect, dxyz)
    ccall((:ocsmGetTessVel, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Ptr{Cdouble}}), modl, ibody, seltype, iselect, dxyz)
end

function ocsmAdjoint(modl, ibody, ndp, ipmtr, irow, icol, nobj, dOdX, dOdD)
    ccall((:ocsmAdjoint, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), modl, ibody, ndp, ipmtr, irow, icol, nobj, dOdX, dOdD)
end

function ocsmTraceStors(modl, pattern, info)
    ccall((:ocsmTraceStors, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), modl, pattern, info)
end

function ocsmGetBody(modl, ibody, type, ichld, ileft, irite, vals, nnode, nedge, nface)
    ccall((:ocsmGetBody, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), modl, ibody, type, ichld, ileft, irite, vals, nnode, nedge, nface)
end

function ocsmPrintBodys(modl, filename)
    ccall((:ocsmPrintBodys, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}), modl, filename)
end

function ocsmPrintBrep(modl, ibody, filename)
    ccall((:ocsmPrintBrep, "libocsm"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cchar}), modl, ibody, filename)
end

function ocsmEvalExpr(modl, expr, value, dot, str)
    ccall((:ocsmEvalExpr, "libocsm"), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}), modl, expr, value, dot, str)
end

function ocsmPrintEgo(obj)
    ccall((:ocsmPrintEgo, "libocsm"), Cvoid, (ego,), obj)
end

function ocsmGetText(icode)
    ccall((:ocsmGetText, "libocsm"), Ptr{Cchar}, (Cint,), icode)
end

function ocsmGetCode(text)
    ccall((:ocsmGetCode, "libocsm"), Cint, (Ptr{Cchar},), text)
end

struct var"union (unnamed at EngSketchPad/include/capsTypes.h:323:3)"
    data::NTuple{56, UInt8}
end

function Base.getproperty(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:323:3)"}, f::Symbol)
    f === :integer && return Ptr{Cint}(x + 0)
    f === :real && return Ptr{Cdouble}(x + 0)
    f === :string && return Ptr{Ptr{Cchar}}(x + 0)
    f === :strings && return Ptr{Ptr{Ptr{Cchar}}}(x + 0)
    f === :tuple && return Ptr{Ptr{capsTuple}}(x + 0)
    f === :pointer && return Ptr{Ptr{Cvoid}}(x + 0)
    f === :own && return Ptr{capsOwn}(x + 0)
    f === :owns && return Ptr{Ptr{capsOwn}}(x + 0)
    f === :errs && return Ptr{Ptr{capsErrs}}(x + 0)
    f === :obj && return Ptr{Ptr{capsObject}}(x + 0)
    f === :objs && return Ptr{Ptr{Ptr{capsObject}}}(x + 0)
    f === :model && return Ptr{ego}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"union (unnamed at EngSketchPad/include/capsTypes.h:323:3)", f::Symbol)
    r = Ref{var"union (unnamed at EngSketchPad/include/capsTypes.h:323:3)"}(x)
    ptr = Base.unsafe_convert(Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:323:3)"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"union (unnamed at EngSketchPad/include/capsTypes.h:323:3)"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

const EGADS_EXTRAPOL = -37

const EGADS_EFFCTOBJ = -36

const EGADS_UVMAP = -35

const EGADS_SEQUERR = -34

const EGADS_CNTXTHRD = -33

const EGADS_READERR = -32

const EGADS_TESSTATE = -31

const EGADS_EXISTS = -30

const EGADS_ATTRERR = -29

const EGADS_TOPOCNT = -28

const EGADS_OCSEGFLT = -27

const EGADS_BADSCALE = -26

const EGADS_NOTORTHO = -25

const EGADS_DEGEN = -24

const EGADS_CONSTERR = -23

const EGADS_TOPOERR = -22

const EGADS_GEOMERR = -21

const EGADS_NOTBODY = -20

const EGADS_WRITERR = -19

const EGADS_NOTMODEL = -18

const EGADS_NOLOAD = -17

const EGADS_RANGERR = -16

const EGADS_NOTGEOM = -15

const EGADS_NOTTESS = -14

const EGADS_EMPTY = -13

const EGADS_NOTTOPO = -12

const EGADS_REFERCE = -11

const EGADS_NOTXFORM = -10

const EGADS_NOTCNTX = -9

const EGADS_MIXCNTX = -8

const EGADS_NODATA = -7

const EGADS_NONAME = -6

const EGADS_INDEXERR = -5

const EGADS_MALLOC = -4

const EGADS_NOTOBJ = -3

const EGADS_NULLOBJ = -2

const EGADS_NOTFOUND = -1

const EGADS_SUCCESS = 0

const EGADS_OUTSIDE = 1

const EGADSMAJOR = 1

const EGADSMINOR = 23

# const EGADSPROP = EGADSprop:Revision(1.23)

const MAGIC = 98789

const MTESSPARAM = 2

const CONTXT = 0

const TRANSFORM = 1

const TESSELLATION = 2

const NIL = 3

const EMPTY = 4

const REFERENCE = 5

const PCURVE = 10

const CURVE = 11

const SURFACE = 12

const NODE = 20

const EDGE = 21

const LOOP = 22

const FACE = 23

const SHELL = 24

const BODY = 25

const MODEL = 26

const EEDGE = 31

const ELOOPX = 32

const EFACE = 33

const ESHELL = 34

const EBODY = 35

const LINE = 1

const CIRCLE = 2

const ELLIPSE = 3

const PARABOLA = 4

const HYPERBOLA = 5

const TRIMMED = 6

const BEZIER = 7

const BSPLINE = 8

const OFFSET = 9

const PLANE = 1

const SPHERICAL = 2

const CYLINDRICAL = 3

const REVOLUTION = 4

const TOROIDAL = 5

const CONICAL = 10

const EXTRUSION = 11

const SREVERSE = -1

const SINNER = -1

const NOMTYPE = 0

const SFORWARD = 1

const SOUTER = 1

const ONENODE = 1

const TWONODE = 2

const OPEN = 3

const CLOSED = 4

const DEGENERATE = 5

const WIREBODY = 6

const FACEBODY = 7

const SHEETBODY = 8

const SOLIDBODY = 9

const ATTRINT = 1

const ATTRREAL = 2

const ATTRSTRING = 3

const ATTRCSYS = 12

const ATTRPTR = 13

const SUBTRACTION = 1

const INTERSECTION = 2

const FUSION = 3

const SPLITTER = 4

const BOX = 1

const SPHERE = 2

const CONE = 3

const CYLINDER = 4

const TORUS = 5

const UISO = 0

const VISO = 1

const NODEOFF = 1

const EDGEOFF = 2

const FACEDUP = 3

const FACECUT = 4

const FACEOFF = 5

const CAPS_SUCCESS = 0

const CAPS_BADRANK = -301

const CAPS_BADDSETNAME = -302

const CAPS_NOTFOUND = -303

const CAPS_BADINDEX = -304

const CAPS_NOTCHANGED = -305

const CAPS_BADTYPE = -306

const CAPS_NULLVALUE = -307

const CAPS_NULLNAME = -308

const CAPS_NULLOBJ = -309

const CAPS_BADOBJECT = -310

const CAPS_BADVALUE = -311

const CAPS_PARAMBNDERR = -312

const CAPS_NOTCONNECT = -313

const CAPS_NOTPARMTRIC = -314

const CAPS_READONLYERR = -315

const CAPS_FIXEDLEN = -316

const CAPS_BADNAME = -317

const CAPS_BADMETHOD = -318

const CAPS_CIRCULARLINK = -319

const CAPS_UNITERR = -320

const CAPS_NULLBLIND = -321

const CAPS_SHAPEERR = -322

const CAPS_LINKERR = -323

const CAPS_MISMATCH = -324

const CAPS_NOTPROBLEM = -325

const CAPS_RANGEERR = -326

const CAPS_DIRTY = -327

const CAPS_HIERARCHERR = -328

const CAPS_STATEERR = -329

const CAPS_SOURCEERR = -330

const CAPS_EXISTS = -331

const CAPS_IOERR = -332

const CAPS_DIRERR = -333

const CAPS_NOTIMPLEMENT = -334

const CAPS_EXECERR = -335

const CAPS_CLEAN = -336

const CAPS_BADINTENT = -337

const CAPS_BADINIT = -338

const CAPS_NOTNEEDED = -339

const CAPS_NOSENSITVTY = -340

const CAPS_NOBODIES = -341

const CAPS_JOURNAL = -342

const CAPS_JOURNALERR = -343

const CAPS_FILELINKERR = -344

# Skipping MacroDefinition: __ProtoExt__ extern

const CAPSLONG = Culong

const DLL = Ptr{Cvoid}

const CAPSMAJOR = 1

const CAPSMINOR = 23

# const CAPSPROP = CAPSprop:Revision(1.23)

const CAPSMAGIC = 1234321

const MAXANAL = 64

const MAXWRITER = 16

const AIM_UPDATESTATE = 1

const AIM_PREANALYSIS = 2

const AIM_POSTANALYSIS = 3

const OCSM_MAJOR_VERSION = 1

const OCSM_MINOR_VERSION = 23

const MAX_NAME_LEN = 64

const MAX_EXPR_LEN = 512

const MAX_FILENAME_LEN = 256

const MAX_LINE_LEN = 2048

const MAX_STR_LEN = 4096

const MAX_STRVAL_LEN = 256

const MAX_STACK_SIZE = 4096

const MAX_NESTING = 20

const MAX_SKETCH_SIZE = 1024

const MAX_SOLVER_SIZE = 256

const MAX_NUM_SKETCHES = 100

const MAX_NUM_MACROS = 100

const SUCCESS = 0

const BAD_MALLOC = -900

const PI = 3.141592653589793

const TWOPI = 6.283185307179586

const PIo2 = 1.5707963267948966

const PIo4 = 0.7853981633974483

const PIo180 = 0.017453292519943295

const HUGEQ = 9.9999999e7

const HUGEI = 9999999

const EPS03 = 0.001

const EPS06 = 1.0e-6

const EPS09 = 1.0e-9

const EPS12 = 1.0e-12

const EPS20 = 1.0e-20

const OCSM_DIMENSION = 100

const OCSM_LBOUND = 101

const OCSM_UBOUND = 102

const OCSM_NAME = 103

const OCSM_ATTRIBUTE = 104

const OCSM_CSYSTEM = 105

const OCSM_POINT = 111

const OCSM_BOX = 112

const OCSM_SPHERE = 113

const OCSM_CONE = 114

const OCSM_CYLINDER = 115

const OCSM_TORUS = 116

const OCSM_IMPORT = 117

const OCSM_UDPRIM = 118

const OCSM_RESTORE = 119

const OCSM_EXTRUDE = 121

const OCSM_RULE = 122

const OCSM_LOFT = 123

const OCSM_BLEND = 124

const OCSM_REVOLVE = 125

const OCSM_SWEEP = 126

const OCSM_FILLET = 131

const OCSM_CHAMFER = 132

const OCSM_HOLLOW = 133

const OCSM_CONNECT = 134

const OCSM_INTERSECT = 141

const OCSM_SUBTRACT = 142

const OCSM_UNION = 143

const OCSM_JOIN = 144

const OCSM_EXTRACT = 145

const OCSM_ELEVATE = 146

const OCSM_TRANSLATE = 151

const OCSM_ROTATEX = 152

const OCSM_ROTATEY = 153

const OCSM_ROTATEZ = 154

const OCSM_SCALE = 155

const OCSM_MIRROR = 156

const OCSM_APPLYCSYS = 157

const OCSM_REORDER = 158

const OCSM_SKBEG = 160

const OCSM_SKVAR = 161

const OCSM_SKCON = 162

const OCSM_LINSEG = 163

const OCSM_CIRARC = 164

const OCSM_ARC = 165

const OCSM_ELLARC = 166

const OCSM_SPLINE = 167

const OCSM_SSLOPE = 168

const OCSM_BEZIER = 169

const OCSM_SKEND = 170

const OCSM_SOLBEG = 171

const OCSM_SOLCON = 172

const OCSM_SOLEND = 173

const OCSM_INTERFACE = 174

const OCSM_END = 175

const OCSM_SET = 176

const OCSM_EVALUATE = 177

const OCSM_GETATTR = 178

const OCSM_UDPARG = 179

const OCSM_SELECT = 180

const OCSM_PROJECT = 181

const OCSM_MACBEG = 182

const OCSM_MACEND = 183

const OCSM_RECALL = 184

const OCSM_STORE = 185

const OCSM_PATBEG = 186

const OCSM_PATBREAK = 187

const OCSM_PATEND = 188

const OCSM_IFTHEN = 189

const OCSM_ELSEIF = 190

const OCSM_ELSE = 191

const OCSM_ENDIF = 192

const OCSM_THROW = 193

const OCSM_CATBEG = 194

const OCSM_CATEND = 195

const OCSM_MARK = 196

const OCSM_GROUP = 197

const OCSM_DUMP = 198

const OCSM_ASSERT = 199

const OCSM_MESSAGE = 200

const OCSM_SPECIAL = 299

const OCSM_PRIMITIVE = 201

const OCSM_GROWN = 202

const OCSM_APPLIED = 203

const OCSM_BOOLEAN = 204

const OCSM_TRANSFORM = 205

const OCSM_SKETCH = 206

const OCSM_SOLVER = 207

const OCSM_UTILITY = 208

const OCSM_ACTIVE = 300

const OCSM_SUPPRESSED = 301

const OCSM_INACTIVE = 302

const OCSM_DEFERRED = 303

const OCSM_SOLID_BODY = 400

const OCSM_SHEET_BODY = 401

const OCSM_WIRE_BODY = 402

const OCSM_NODE_BODY = 403

const OCSM_NULL_BODY = 404

const OCSM_DESPMTR = 500

const OCSM_CFGPMTR = 501

const OCSM_CONPMTR = 502

const OCSM_LOCALVAR = 503

const OCSM_OUTPMTR = 504

const OCSM_UNKNOWN = 505

const OCSM_NODE = 600

const OCSM_EDGE = 601

const OCSM_FACE = 602

const OCSM_BODY = 603

const OCSM_EEDGE = 604

const OCSM_EFACE = 605

const OCSM_EBODY = 606

const OCSM_UNDEFINED = -123.456

const OCSM_FILE_NOT_FOUND = -201

const OCSM_ILLEGAL_STATEMENT = -202

const OCSM_NOT_ENOUGH_ARGS = -203

const OCSM_NAME_ALREADY_DEFINED = -204

const OCSM_NESTED_TOO_DEEPLY = -205

const OCSM_IMPROPER_NESTING = -206

const OCSM_NESTING_NOT_CLOSED = -207

const OCSM_NOT_MODL_STRUCTURE = -208

const OCSM_PROBLEM_CREATING_PERTURB = -209

const OCSM_EBODY_NOT_FOUND = -210

const OCSM_MISSING_MARK = -211

const OCSM_INSUFFICIENT_BODYS_ON_STACK = -212

const OCSM_WRONG_TYPES_ON_STACK = -213

const OCSM_DID_NOT_CREATE_BODY = -214

const OCSM_CREATED_TOO_MANY_BODYS = -215

const OCSM_TOO_MANY_BODYS_ON_STACK = -216

const OCSM_ERROR_IN_BODYS_ON_STACK = -217

const OCSM_MODL_NOT_CHECKED = -218

const OCSM_NEED_TESSELLATION = -219

const OCSM_BODY_NOT_FOUND = -221

const OCSM_FACE_NOT_FOUND = -222

const OCSM_EDGE_NOT_FOUND = -223

const OCSM_NODE_NOT_FOUND = -224

const OCSM_ILLEGAL_VALUE = -225

const OCSM_ILLEGAL_ATTRIBUTE = -226

const OCSM_ILLEGAL_CSYSTEM = -227

const OCSM_NO_SELECTION = -228

const OCSM_SKETCH_IS_OPEN = -231

const OCSM_SKETCH_IS_NOT_OPEN = -232

const OCSM_COLINEAR_SKETCH_POINTS = -233

const OCSM_NON_COPLANAR_SKETCH_POINTS = -234

const OCSM_TOO_MANY_SKETCH_POINTS = -235

const OCSM_TOO_FEW_SPLINE_POINTS = -236

const OCSM_SKETCH_DOES_NOT_CLOSE = -237

const OCSM_SELF_INTERSECTING = -238

const OCSM_ASSERT_FAILED = -239

const OCSM_ILLEGAL_CHAR_IN_EXPR = -241

const OCSM_CLOSE_BEFORE_OPEN = -242

const OCSM_MISSING_CLOSE = -243

const OCSM_ILLEGAL_TOKEN_SEQUENCE = -244

const OCSM_ILLEGAL_NUMBER = -245

const OCSM_ILLEGAL_PMTR_NAME = -246

const OCSM_ILLEGAL_FUNC_NAME = -247

const OCSM_ILLEGAL_TYPE = -248

const OCSM_ILLEGAL_NARG = -249

const OCSM_NAME_NOT_FOUND = -251

const OCSM_NAME_NOT_UNIQUE = -252

const OCSM_PMTR_IS_DESPMTR = -253

const OCSM_PMTR_IS_LOCALVAR = -254

const OCSM_PMTR_IS_OUTPMTR = -255

const OCSM_PMTR_IS_CONPMTR = -256

const OCSM_WRONG_PMTR_TYPE = -257

const OCSM_FUNC_ARG_OUT_OF_BOUNDS = -258

const OCSM_VAL_STACK_UNDERFLOW = -259

const OCSM_VAL_STACK_OVERFLOW = -260

const OCSM_ILLEGAL_BRCH_INDEX = -261

const OCSM_ILLEGAL_PMTR_INDEX = -262

const OCSM_ILLEGAL_BODY_INDEX = -263

const OCSM_ILLEGAL_ARG_INDEX = -264

const OCSM_ILLEGAL_ACTIVITY = -265

const OCSM_ILLEGAL_MACRO_INDEX = -266

const OCSM_ILLEGAL_ARGUMENT = -267

const OCSM_CANNOT_BE_SUPPRESSED = -268

const OCSM_STORAGE_ALREADY_USED = -269

const OCSM_NOTHING_PREVIOUSLY_STORED = -270

const OCSM_SOLVER_IS_OPEN = -271

const OCSM_SOLVER_IS_NOT_OPEN = -272

const OCSM_TOO_MANY_SOLVER_VARS = -273

const OCSM_UNDERCONSTRAINED = -274

const OCSM_OVERCONSTRAINED = -275

const OCSM_SINGULAR_MATRIX = -276

const OCSM_NOT_CONVERGED = -277

const OCSM_UDP_ERROR1 = -281

const OCSM_UDP_ERROR2 = -282

const OCSM_UDP_ERROR3 = -283

const OCSM_UDP_ERROR4 = -284

const OCSM_UDP_ERROR5 = -285

const OCSM_UDP_ERROR6 = -286

const OCSM_UDP_ERROR7 = -287

const OCSM_UDP_ERROR8 = -288

const OCSM_UDP_ERROR9 = -289

const OCSM_OP_STACK_UNDERFLOW = -291

const OCSM_OP_STACK_OVERFLOW = -292

const OCSM_RPN_STACK_UNDERFLOW = -293

const OCSM_RPN_STACK_OVERFLOW = -294

const OCSM_TOKEN_STACK_UNDERFLOW = -295

const OCSM_TOKEN_STACK_OVERFLOW = -296

const OCSM_UNSUPPORTED = -298

const OCSM_INTERNAL_ERROR = -299

const OCSM_EXTERNAL = OCSM_DESPMTR

const OCSM_CONFIG = OCSM_CFGPMTR

const OCSM_CONSTANT = OCSM_CONPMTR

const OCSM_OUTPUT = OCSM_OUTPMTR

const OCSM_INTERNAL = OCSM_LOCALVAR

const OCSM_PMTR_IS_EXTERNAL = OCSM_PMTR_IS_DESPMTR

const OCSM_PMTR_IS_CONSTANT = OCSM_PMTR_IS_CONPMTR

const OCSM_PMTR_IS_OUTPUT = OCSM_PMTR_IS_OUTPMTR

const OCSM_PMTR_IS_INTERNAL = OCSM_PMTR_IS_LOCALVAR

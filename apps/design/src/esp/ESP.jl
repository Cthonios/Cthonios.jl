module ESP

export CoordinateSensitivity
export Geometry
export Parameters
export Tesselation
export tesselate!

using CEnum
using EngineeringSketchPad_jll
using Exodus
using Exodus_jll
using OCCT_jll
using Printf

struct verTags
    ptype::Cint
    pindex::Cint
end

function _check_egads_status(status, msg=nothing)
    if status != 0
        if msg !== nothing
            _egads_err(status)
        end
    end
end

function _check_ocsm_status(status, msg=nothing)
    if status != 0
        if msg !== nothing
            @error msg
        end
        error("EGADS ERROR CODE: $status")
    end
end

@inline _devnull_path() = Sys.iswindows() ? "NUL" : "/dev/null"

@inline function _egads_err(status)
    print(stderr, "EGADS ERROR CODE: ")
    print(stderr, status)
    error()
end

include("Wrappers.jl")

include("Parameters.jl")
include("Geometry.jl")
include("Tesselation.jl")
include("Sensitivity.jl")

function __init__()
    esp_jll_dir = EngineeringSketchPad_jll.artifact_dir
    ENV["ESP_ARCH"] = "LINUX64"
    ENV["ESP_ROOT"] = esp_jll_dir
    ENV["CASROOT"] = OCCT_jll.artifact_dir
    ENV["CARARCH"] = "."
    ENV["CASREV"] = "7.7"
    ENV["CAPS_GLYPH"] = joinpath(esp_jll_dir, "src", "CAPS", "aim", "pointwise", "glyph")
    ENV["ESP_JULIA_ROOT"] = dirname(dirname(@__FILE__))
    # collect all dep libraries
    dep_libs = Exodus_jll.LIBPATH[] * EngineeringSketchPad_jll.LIBPATH[]
    ENV["PYESPDEPLIBS"] = dep_libs
    ENV["UDUNITS2_XML_PATH"] = joinpath(esp_jll_dir, "share", "esp_udunits", "udunits2.xml")
end

end # module ESP

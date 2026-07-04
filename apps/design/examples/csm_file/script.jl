import CthoniosDesign.ESP as ESP
using CthoniosDesign.ESP

csm_file = Base.source_dir() * "/hole.csm"
geo = Geometry(csm_file; mode = "julia", verbose = false)
params = copy(get_design_parameters(geo))
params[1] = 2.5
ESP.update!(geo, params; verbose = false)
# ESP.dump_sensitivities(geo, "sens.exo")

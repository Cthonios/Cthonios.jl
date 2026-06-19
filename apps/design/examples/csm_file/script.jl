using CthoniosDesign.ESP

csm_file = Base.source_dir() * "/hole.csm"
sensitivity = CoordinateSensitivity(csm_file)

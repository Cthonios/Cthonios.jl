import CthoniosDesign as CD

csm_file = Base.source_dir() * "/hole.csm"
opt = CD.StructuralOptimization(csm_file; esp_mode = "julia")

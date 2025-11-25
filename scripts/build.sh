julia \
    --project=@. \
    scripts/CthoniosBuild.jl \
    --build-sysimage
julia \
    --project=@. \
    --sysimage cthonios_sysimage.so \
    scripts/CthoniosBuild.jl \
    --build-app \
    --incremental

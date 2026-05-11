# using JuliaC
# build_path = joinpath(dirname(@__DIR__), "build")
# src_path = joinpath(dirname(@__DIR__))
# @show build_path
# @show src_path
# rm(build_path; force=true, recursive=true)

# img = ImageRecipe(
#     output_type    = "--output-exe",
#     file           = "$src_path",
#     trim_mode      = nothing,
#     add_ccallables = false,
#     verbose        = false,
# )

# link = LinkRecipe(
#     image_recipe = img,
#     outname      = "$build_path/cthonios"
# )

# bun = BundleRecipe(
#     link_recipe = link,
#     output_dir  = build_path # or `nothing` to skip bundling
# )

# compile_products(img)
# link_products(link)
# bundle_products(bun)

using PackageCompiler
src_path = joinpath(dirname(@__DIR__))
build_path = joinpath(dirname(@__DIR__), "build")
create_app(src_path, build_path; force = true)

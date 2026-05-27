using JuliaC
build_path = joinpath(dirname(@__DIR__), "build")
src_path = joinpath(dirname(@__DIR__))
@show build_path
@show src_path
rm(build_path; force=true, recursive=true)

img = ImageRecipe(
    output_type    = "--output-exe",
    file           = "$src_path/src/Main.jl",
    trim_mode      = "safe",
    add_ccallables = false,
    verbose        = false,
)

link = LinkRecipe(
    image_recipe = img,
    outname      = "$build_path/cthonios"
)

bun = BundleRecipe(
    link_recipe = link,
    output_dir  = build_path # or `nothing` to skip bundling
)

compile_products(img)
link_products(link)
bundle_products(bun)

using ArgParse
using PackageCompiler
using TOML

@info "Running CthoniosBuild.jl"

settings = ArgParseSettings()
@add_arg_table! settings begin
  "--clean", "-c"
    help = "Clean the build directory"
    action = :store_true
  "--build-dir"
    help = "Directory to build sysimage and executables"
    arg_type = String
    default = "./build"
  "--build-sysimage"
    help = "Build a sysimage for the current project"
    action = :store_true
  "--build-app"
    help = "Build an app in build-dir"
    action = :store_true
  "--incremental"
    help = "Use an incremental sysimage"
    action = :store_true
end

settings = parse_args(settings)

if settings["clean"]
  @info "Cleaning build directory"
  rm(settings["build-dir"]; force=true, recursive=true)
  exit()
end

if settings["build-sysimage"]
  @info "Building sysimage"
  @info "Parsing deps from Project.toml"
  @info "Found the following deps"
  toml = TOML.parsefile("Project.toml")
  deps = keys(toml["deps"]) |> collect
  deps = sort(deps)
  display(deps)
  # TODO figure out how to use weakdeps
  # @info "Found the following weakdeps"
  # weakdeps = keys(toml["weakdeps"])
  # display(weakdeps)
  # all_deps = vcat(deps..., weakdeps...) |> sort
  # @info "Using all the following deps and weakdeps"
  # display(all_deps)
  create_sysimage(deps, sysimage_path="cthonios_sysimage.so"; cpu_target="generic")
end

if settings["build-app"]
  @info "Building app"
  if settings["incremental"]
    @info "Using incremental sysimage"
    
    incremental = true
  else
    incremental = false
  end
  create_app(
    "./", settings["build-dir"],
    executables=executables=["cthonios" => "cthonios_main"],
    force=true;
    cpu_target="generic",
    # base_sysimage="cthonios_sysimage.so"
    include_lazy_artifacts=true,
    incremental=incremental
  )
end

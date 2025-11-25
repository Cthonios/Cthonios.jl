using ArgParse
using PackageCompiler
using Pkg
using TOML

@info "Running CthoniosBuild.jl"

settings = ArgParseSettings()
@add_arg_table! settings begin
  "--build-app"
    help = "Build an app in build-dir"
    action = :store_true
  "--build-dir"
    help = "Directory to build sysimage and executables"
    arg_type = String
    default = "./build"
  "--build-sysimage"
    help = "Build a sysimage for the current project"
    action = :store_true
  "--clean", "-c"
    help = "Clean the build directory"
    action = :store_true
  "--include-all-extensions"
    help = "Build all available extensions"
    action = :store_true
  "--incremental"
    help = "Use an incremental sysimage"
    action = :store_true
end

settings = parse_args(settings)

@info "Checking project status - note you should use a temporary environment"
Pkg.status()

if settings["clean"]
  @info "Cleaning build directory"
  rm(settings["build-dir"]; force=true, recursive=true)
  exit()
end

if settings["build-sysimage"]
  @info "Building sysimage"
  # @info "Parsing deps from Project.toml"
  # @info "Found the following deps"
  # toml = TOML.parsefile("Project.toml")
  # deps = keys(toml["deps"]) |> collect

  if settings["include-all-extensions"]
    extra_deps = [
      "Adapt",
      "ArgParse",
      "AMDGPU",
      "EngineeringSketchPadWrapper",
      "YAML"
    ]
  else
    # just build the CLI extension
    extra_deps = [
      # "ArgParse",
      # "YAML"
    ]
  end
  
  # deps = sort(deps)
  # display(deps)
  @info "Adding the following weak deps"
  display(extra_deps)
  Pkg.add.(extra_deps)

  # @info "Adding Cthonios.jl as dep"
  # Pkg.develop(path=".")
  # @info "Instantiating fresh temporary environment"
  # Pkg.instantiate()

  # TODO figure out how to use weakdeps
  # @info "Found the following weakdeps"
  # weakdeps = keys(toml["weakdeps"])
  # display(weakdeps)
  # all_deps = vcat(deps..., weakdeps...) |> sort
  # @info "Using all the following deps and weakdeps"
  # display(all_deps)

  @info "Parsing deps from Project.toml"
  @info "Found the following deps"
  toml = TOML.parsefile(dirname(Base.active_project()) * "/Project.toml")
  deps = keys(toml["deps"]) |> collect

  create_sysimage(
    deps, sysimage_path="cthonios_sysimage.so"; 
    cpu_target="generic",
    include_transitive_dependencies=true
  )
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
    incremental=incremental,
    # precompile_execution_file="$(dirname(@__FILE__))/precompile.jl"
  )
end

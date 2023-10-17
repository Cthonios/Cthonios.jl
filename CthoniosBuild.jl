import REPL
using Logging
using PackageCompiler
using REPL.TerminalMenus

options = ["Continue", "Overwrite"]
menu = RadioMenu(options, pagesize=2)

function check_for_build()
  # println(">>> Checking for existing Cthonios build...")
  @info "Checking for existing Cthonios build..."
  build_dir = abspath(ENV["CTHONIOS_BUILD_DIR"])
  if isdir(build_dir)
    choice = request("Existing build found in $build_dir. Overwrite build?", menu)

    if choice != -1
      if choice == 1
        # println(">>> Continuing CthoniosBuild.jl with existing buld directory...")
        @info "Continuing CthoniosBuild.jl with existing buld directory..."
        # exit()
        return
      elseif choice == 2
        # println(">>> Removing old build in $build_dir...")
      @info "Removing old build in $build_dir..."

        rm(build_dir; force=true, recursive=true)
      end
    else
      # println("Menu cancelled")
      @info "Menu cancelled"
    end
  end

  # println(">>> Making new build directory in $build_dir...")
  @info "Making new build directory in $build_dir..."
  mkdir(build_dir)
end

function build_sysimage()
  # println(">>> Checking for existing sysimage...")
  @info "Checking for existing sysimage..."
  build_dir = ENV["CTHONIOS_BUILD_DIR"]
  sysimage_path = joinpath(build_dir, "sysimage.so")
  if isfile(sysimage_path)
    choice = request("Existing sysimage found in $build_dir. Overwrite sysimage?", menu)

    if choice != -1
      if choice == 1
        # println(">>> Continuing CthoniosBuild.jl with existing sysimage...")
        @info "Continuing CthoniosBuild.jl with existing sysimage..."
        # exit()
        return
      elseif choice == 2
        # println(">>> Removing old sysimage...")
        @info "Removing old sysimage..."
        rm(sysimage_path; force=true)
      end
    else
      println("Menu cancelled")
    end  
  end

  # println(">>> Building sysimage in $(abspath(build_dir))...")
  @info "Building sysimage in $(abspath(build_dir))..."
  # TODO just read in deps from Project.toml
  create_sysimage(
    [
      "ArgParse",
      "Exodus",
      "FiniteElementContainers",
      "IterativeSolvers",
      "Logging",
      "LoggingExtras",
      "ReferenceFiniteElements",
      "YAML"
    ];
    sysimage_path=sysimage_path,
    cpu_target="x86-64",
    incremental=true
  )
end

function build_executable()
  # println(">>> Checking for existing Cthonios executable build dir...")
  @info "Checking for existing Cthonios executable build dir..."
  build_dir = joinpath(ENV["CTHONIOS_BUILD_DIR"], "build_dir")

  create_app(
    "./", build_dir;
    executables=["cthonios" => "julia_main"],
    incremental=true
  )
end

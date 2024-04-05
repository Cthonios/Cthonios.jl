using PackageCompiler
create_app(
  "./", "cthonios";
  executables=[
    "cthonios" => "julia_main"
  ],
  # filter_stdlibs=true # segfault barf
  # precompile_execution_file="precompile/precompile.jl"
)

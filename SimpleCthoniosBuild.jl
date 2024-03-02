using PackageCompiler
create_app(
  "./", "cthonios";
  executables=[
    "cthonios" => "julia_main"
  ],
  precompile_execution_file="precompile/precompile.jl"
)

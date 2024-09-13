using PackageCompiler
create_app(
  "./", "build"; 
  executables=executables=["cthonios" => "cthonios_main"],
  force=true
)

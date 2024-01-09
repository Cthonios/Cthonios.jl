CPU_NAME='x86_64' julia --project=@. -e '
  using PackageCompiler
  create_app("./", "cthonios")
'

using Cthonios

for dir in readdir("precompile", join=true)
  if isdir(dir)
    args = [
      "-i", joinpath(dir, "precompile.yaml")
    ]
    append!(ARGS, args)
    Cthonios.julia_main()
    pop!(ARGS)
    pop!(ARGS)
    rm(joinpath(dir, "precompile.log"))
    rm(joinpath(dir, "precompile.e"))
  end
end

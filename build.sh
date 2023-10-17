CTHONIOS_BUILD_DIR=./build julia --project=@. -e '
  include("CthoniosBuild.jl")
  check_for_build()
  build_sysimage()
'
CTHONIOS_BUILD_DIR=./build julia --project=@. --sysimage=./build/sysimage.so -e '
  include("CthoniosBuild.jl")
  build_executable()
'

```@meta
CurrentModule = Cthonios
```
# Running from the REPL
To run ```Cthonios``` from the REPL one can proceed in one of two manners. Either you can write your input to ```Cthonios``` as a julia script (e.g. ```script.jl```) or you can provide input via our ```yaml``` input file.

The script approach can proceed with the following recipe
```julia
using Revise, Cthonios
include("script.jl")
```

The input file approach can proceed with the following recipe
```julia
using Revise, Cthonios
push!(ARGS, "-i")
push!(ARGS, "input_file.yaml")
Cthonios.cthonios_main()
```

Of course one can just run the following from a terminal
```
julia script.jl
```
but this will re-compile certain features of the code for each run which is not great for de-bugging a new problem.

# Running Cthonios as an executable
To run ```Cthonios``` as an executable, one must first run the ```CthoniosBuild.jl``` script (assuming one has ```PackageCompiler``` installed already). This will build an executable called ```cthonios``` in a ```build``` folder. This script has several options and we suggest running 
```
julia --project=@. CthoniosBuild.jl -h
```
prior to proceeding with an executable.

To run ```Cthonios``` as an executable you can use the following command (after building of course)

```
/path/to/cthonios -i <input-file.yaml>
```

See the examples folder for input file examples. Please note that input file syntax will usually lag behind newly implemented capabilites since we do not strictly enforce each feature to implement yaml parsing currently.

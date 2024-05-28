# Different function types
"""
$(TYPEDEF)
"""
const ScalarFunction{D, T1, T2, T3} = FunctionWrapper{T1, Tuple{SVector{D, T2}, T3}}

# Different containers for function types
const ScalarFunctionWrapperArray = Vector{ScalarFunction}

"""
Method to define a scalar function via a script
"""
function ScalarFunction(f::F, n_dims::Int) where F <: Function
  return ScalarFunction{n_dims, Float64, Float64, Float64}(f)
end

"""
$(TYPEDSIGNATURES)
This method constructs a ```Function``` from a
a dictionary of inputs. Typically this will be called
from ```setup_functions``` which will recieve a
vector of dictionaries from a parsed yaml file.

The syntax in the input file should mimic the following
```
bc_name:
  type: ScalarFunction{NDIM, Float64, Float64, Float64}
  expression: (x, t) -> 1.0 * t
```
"""
function setup_function(inputs::D) where D <: Dict{Symbol, Any}
  type = inputs[:type]
  func = inputs[:expression]
  @info "Function"
  @info "  type       = $type"
  @info "  expression = $func"
  @info ""
  # return eval(Meta.parse(type))(eval(Meta.parse(func)))
  @warn "Currently not using ScalarFunction..."
  ex = Meta.parse(func)
  return @RuntimeGeneratedFunction(ex)
end

"""
$(TYPEDSIGNATURES)
Input file syntax

```
functions:
  bc_name_1:
    type: ScalarFunction{NDIM, Float64, Float64, Float64}
    expression: (x, t) -> 1.0 * t
  bc_name_2:
    type: ScalarFunction{NDIM, Float64, Float64, Float64}
    expression: (x, t) -> 1.0 * x * t
  ...
  bc_name_n:
    type: ScalarFunction{NDIM, Float64, Float64, Float64}
    expression: (x, t) -> 0.0
```
"""
function setup_functions(inputs::D) where D <: Vector{Dict{Symbol, Any}}
  new_section("Functions")
  funcs = map(x -> setup_function(x), inputs)
  return funcs
end

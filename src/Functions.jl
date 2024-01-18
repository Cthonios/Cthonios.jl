# Different function types
const ScalarFunction{D, T1, T2, T3} = FunctionWrapper{T1, Tuple{SVector{D, T2}, T3}}

# Different containers for function types
const ScalarFunctionWrapperArray = Vector{ScalarFunction}

function setup_function(inputs::D) where D <: Dict{Symbol, Any}
  type = inputs[:type]
  func = inputs[:expression]
  @info "Function"
  @info "  type       = $type"
  @info "  expression = $func"
  @info ""
  return eval(Meta.parse(type))(eval(Meta.parse(func)))
end

function setup_functions(inputs::D) where D <: Vector{Dict{Symbol, Any}}
  new_section("Functions")
  funcs = map(x -> setup_function(x), inputs)
  return funcs
end

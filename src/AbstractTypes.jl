abstract type AbstractCthoniosType end

function from_yaml(type::AbstractCthoniosType, dict)
  @assert "need to overload this if you want to parse"
end

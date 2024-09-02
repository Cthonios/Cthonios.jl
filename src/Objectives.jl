"""
$(TYPEDEF)
"""
abstract type AbstractObjective end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Objective{D, F1, F2, F3}
  domain::D
  value::F1
  gradient::F2
  hessian::F3
end

# value(o::Objective, Uu, p) = o.value(o.domain, Uu, p)
"""
$(TYPEDSIGNATURES)
"""
gradient(o::Objective, U, p) = o.gradient(o.domain, U, p)
"""
$(TYPEDSIGNATURES)
"""
hessian(o::Objective, U, p) = o.hessian(o.domain, U, p)


# exports
export Objective
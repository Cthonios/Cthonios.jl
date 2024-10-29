struct Source{NF, F} <: AbstractPhysics{NF, 0, 0}
  func::F
end

function Source{NF}(func) where NF
  return Source{NF, typeof(func)}(func)
end

function energy(physics::Source, u::T, ∇u, X, t, Z, props) where T
  b = physics.func(X, current_time(t))
  return dot(b, u)
end

function gradient(physics::Source, u, ∇u, v, ∇v, X, t, Z, props)
  b = physics.func(X, current_time(t))
  R = b * v'
  return R
end

function hessian(physics::Source, u, ∇u, v, ∇v, X, t, Z, props)
  N = length(v) * num_fields(physics)
  return zeros(SMatrix{N, N, eltype(u), N * N})
end

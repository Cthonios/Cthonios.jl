struct Dynamics{NF} <: AbstractPhysics{NF, 0, 0}
end

function energy(physics::Dynamics, u::T, ∇u, X, t, Z, props) where T
  density = 1.0 # TODO
  return 0.5 * density * dot(u, u)
end

function gradient(physics::Dynamics, u, ∇u, v, ∇v, X, t, Z, props)
  density = 1.0 # TODO
  return density * u * v'
end

function hessian(physics::Dynamics, u, ∇u, v, ∇v, X, t, Z, props)
  density = 1.0 # TODO
  return density * v * v'
end

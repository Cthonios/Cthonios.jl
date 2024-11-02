struct Laplacian{NF} <: AbstractPhysics{NF, 0, 0}
end

function energy(physics::Laplacian, u::T, ∇u, X, t, Z, props) where T
  return 0.5 * dot(∇u, ∇u)
end

function gradient(physics::Laplacian, u, ∇u, v, ∇v, X, t, Z, props)
  R = ∇u * ∇v'
  return R
end

function hessian(physics::Laplacian, u, ∇u, v, ∇v, X, t, Z, props)
  K = ∇v * ∇v'
  return K
end

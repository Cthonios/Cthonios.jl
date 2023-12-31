import AbstractDifferentiation as AD
using Cthonios
using ForwardDiff
using IterativeSolvers
using LinearAlgebra
using LinearMaps
using LinearOperators
using LinearSolve
using Preconditioners
using SciMLOperators
using YAML

function objective_func(domain, X, Uu)
  Cthonios.energy(domain, X, Uu)
end

function LinearAlgebra.ldiv!(y, A::LinearOperator, x)
  y .= A * x
end

backend = AD.ForwardDiffBackend()
domain  = Cthonios.QuasiStaticDomain(YAML.load_file("test/test_input_file.yaml"))

obj = Cthonios.Objective(domain, backend, objective_func)

Cthonios.step!(obj.domain.time)
Cthonios.step!(obj.domain.time)
Cthonios.step!(obj.domain.time)
Cthonios.step!(obj.domain.time)
Cthonios.update_unknown_ids!(obj.domain)

X  = obj.domain.coords
Uu = Cthonios.create_unknowns(obj.domain)
v   = Cthonios.create_unknowns(domain)
v   .= rand(Float64, size(v))
v_x = Cthonios.create_fields(obj.domain)
v_x .= rand(Float64, size(v_x))
Cthonios.objective(obj, X, Uu)
# Cthonios.grad_x(obj, X, Uu)
Cthonios.grad_u(obj, X, Uu)
# @time Cthonios.hvp_x(obj, X, Uu)(v_x)[1]
# @time Cthonios.hvp_u(obj, X, Uu)(v)[1]

# residual
R = Cthonios.grad_u(obj, X, Uu)
K = Cthonios.stiffness(domain, Uu)
# display(K)
hvp = Cthonios.hvp_u(obj, X, Uu)

f(v, p, t) = hvp(v)[1]
op = FunctionOperator(f, zeros(length(Uu)), zeros(length(Uu)))
# cholesky(op)
# CholeskyFactorization(f)
# cholesky(f)

# m = LinearMap(f, length(Uu), length(Uu); issymmetric=true)
Diagonal(op)
# cholesky(m)
# CholeskyFactorization(m)
# m(v)
# K * v
# @show m(v) ≈ K * v
# function mulprecond!(res, v, α, β)
#   res .= α * f(v)
# end

# @show f(v) ≈ K * v

# m = LinearMap(f, length(Uu))
# m * v
# InverseMap(m) * v
# CholeskyPreconditioner(K, 2)
# prob = LinearProblem(op, Uu)
# solve(prob)
# prob = LinearProblem(K, R)
# sol = solve(prob, KrylovJL_GMRES())
# cg(m, R)
# K \ R
# op = LinearOperator(
#   Float64, length(Uu), length(Uu), true, false, 
#   mulprecond!
# )

# op = opCholesky(K)
# norm(op * v - K \ v) / norm(v)
# cholesky(K)
# display(K)
# display(
# op = Preconditioners.CholeskyPreconditioner(K, 1)

# cg(K, R; Pl=op)
# @show K \ R ≈ cg(K, R; Pl=op)
# @time K \ R
# @time cg(K, R; Pl=op)
# K \ v
# K .- K |> norm

# @show op * v == f(v)
# Preconditioners.CholeskyPreconditioner(op, 2)
# opCholesky(op)

# op \ Uu

# hvp(v)[1]
# 
# m = LinearMap(f, length(Uu); issymmetric=true)
# InverseMap(m)
# m = LinearMap(Cthonios.hvp_u(obj, X, Uu), length(Uu); issymmetric=true)
# p = CholeskyPreconditioner(m, 2)
# # InverseMap(m) * Uu

# R = Cthonios.grad_u(obj, X, Uu)
# # @show typeof(R)
# # @show size(R)
# @show size(m)
# @show size(Uu)
# m(Uu)
# InverseMap(m) * R
# m \ R
# IterativeSolvers.cg(m, R)

# Cthonios.func_map(obj, X, Uu)

# Cthonios.energy(obj.domain, Uu)
# U   = Cthonios.create_fields(domain)
# X   = domain.coords



# Uu  = Cthonios.create_unknowns(domain)


# Cthonios.update_bcs!(U, domain)
# Cthonios.update_unknowns!(U, domain, Uu)

# @show objective(domain, Uu)
# @show Cthonios.energy(domain, Uu)


# ad = AD.gradient(backend, x -> objective(domain, x), Uu)[1]
# an = Cthonios.residual(domain, Uu)
# @show ad ≈ an

# hvp = AD.pushforward_function(backend, x -> Cthonios.residual(domain, x), Uu)

# K = Cthonios.stiffness(domain, Uu)
# Kv = K * v
# @time Hv = hvp(v)[1]

# @show Kv ≈ Hv

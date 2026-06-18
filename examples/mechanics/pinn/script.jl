using ConstitutiveModels
using Cthonios
using FiniteElementContainers
using LinearAlgebra
using Lux
using Optimisers
using Random
using Zygote

mesh_file = Base.source_dir() * "/mesh.g"

function solution(asm, x, t, z)
    display(z)
    L = 1.0
    u_max = 0.05

    u = create_field(asm)
    @. u[1, :] = x[1, :] * t * u_max / L + x[1, :] * (x[1, :] - L) * t * z[1, :] / L / L
    @. u[2, :] = x[1, :] * (x[1, :] - L) * t * z[2, :] / L / L
    return u
end

# function gradient(asm, θs, st, p)
#     x = p.coords
#     t = p.times.time_current
#     z, back = Zygote.pullback(θs) do θ
#         y, _ = Lux.apply(model, x, θ, st)
#         y
#     end
#     u = solution(asm, x, t, z)
#     assemble_vector!(asm, residual!, u, p)
#     # dLdu = residual(asm)
#     dLdu = asm.residual_storage
#     display(u)
#     display(dLdu)
#     dθs = back(dLdu)[1]
#     return dθs
# end

function bc_pullback(dLdu, x, t)

    L = 1.0

    φ = @. x[1,:] * (x[1,:] - L) * t / (L^2)

    return reshape(φ, 1, :) .* dLdu
end

function gradient(asm, θs, st, p)

    x = p.coords
    t = p.times.time_current

    z, back = Zygote.pullback(θs) do θ
        y, _ = Lux.apply(model, x, θ, st)
        y
    end
    println("z min = ", minimum(z))
    println("z max = ", maximum(z))

    u = solution(asm, x, t, z)

    assemble_vector!(asm, residual!, u, p)

    dLdu = asm.residual_storage

    dLdz = bc_pullback(dLdu, x, t)

    return back(dLdz)[1]
end

function loss_and_gradient(asm, p, model, θs, st)

    x = p.coords
    t = p.times.time_current

    z, back = Zygote.pullback(θs) do θ
        y, _ = Lux.apply(model, x, θ, st)
        y
    end
    z = H1Field(z)
    # u = solution(asm, x, t, z)
    u = create_field(asm)
    @. u[1,:] = 0.05 * p.coords[1,:]
    @. u[2,:] = 0.0
    # fill!(u, 0.0)
    # @show minimum(u)
    # @show maximum(u)
    display(u)
    # assemble_vector!(asm, residual!, u, p)
    assemble_scalar!(asm, energy, u, p)

    # @show "Hur"
    loss = sum(asm.scalar_quadrature_storage)

    dLdu = asm.residual_storage

    dLdz = bc_pullback(dLdu, x, t)

    grad = back(dLdz)[1]

    return loss, grad
end

function main_func()

    # rng
    rng = Random.default_rng()

    # setup mesh
    mesh = UnstructuredMesh(mesh_file)
    # x = mesh.nodal_coords
    # x = Matrix(x)

    V = FunctionSpace(mesh, H1Field, Lagrange) 
    physics = SolidMechanics(PlaneStrain(), NeoHookean())
    prop_inputs = Dict{String, Any}(
        "density"       => 1.0,
        "bulk modulus"  => 0.833,
        "shear modulus" => 0.3846
    )
    props = create_properties(physics, prop_inputs)
    u = VectorFunction(V, "displ")
    asm = SparseMatrixAssembler(u; use_condensed = true, use_inplace_methods = true)
    times = TimeStepper(0., 1., 1)

    # bcs
    fixed(_, _) = 0.0
    displace(_, t) = 0.05 * t
    dbcs = DirichletBC[
        DirichletBC("displ_x", fixed; nodeset_name = "nset_2"),
        DirichletBC("displ_y", fixed; nodeset_name = "nset_2"),
        DirichletBC("displ_x", displace; nodeset_name = "nset_4"),
        DirichletBC("displ_y", fixed; nodeset_name = "nset_4"),
    ]

    p = create_parameters(mesh, asm, physics, props; dirichlet_bcs = dbcs, times = times)
    FiniteElementContainers.update_time!(p)
    # neural network setup
    model = Chain(
        Dense(2, 100, tanh),
        Dense(100, 100, tanh),
        Dense(100, 2)
    ) |> f64

    θs, st = Lux.setup(rng, model)
    θs = θs |> f64
    st = st |> f64

    θs.layer_3.weight .= 0
    θs.layer_3.bias .= 0

    opt = Optimisers.Adam(1e-3)
    opt_state = Optimisers.setup(opt, θs)

    nepochs = 5000

    for epoch in 1:nepochs
        @show "Hur"
        loss, grad = loss_and_gradient(asm, p, model, θs, st)
        opt_state, θs = Optimisers.update(opt_state, θs, grad)

        if epoch % 100 == 0
            println("epoch = ", epoch, ", loss = ", loss, ", ||g|| = ", norm(residual(asm)))
            display(asm.residual_storage)
        end
    end
end

main_func()

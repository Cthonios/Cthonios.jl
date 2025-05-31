function FiniteElementContainers.residual(
    physics::AbstractPhysics, interps, u_el, x_el, state_old_q, props_el, t, dt
)
    (; X_q, N, ∇N_X, JxW) = MappedInterpolants(interps, x_el)
    u_el = reshape_element_level_field(physics, u_el)
    
end

abstract type AbstractPhysicsKernel{NF, NP, NS} <: AbstractPhysics{NF, NP, NS} end

include("kernels/Laplacian.jl")

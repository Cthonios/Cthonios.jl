abstract type AbstractSolidMechanics{NF, NP, NS, FM, CM} <: AbstractPhysics{NF, NP, NS} end

# this one could be implemented more intelligently
function characteristic_element_length(
    x_el::SMatrix{1, N, T, N2},
    u_el::SMatrix{1, N, T, N2} = zero(u_el::SMatrix{1, N, T, N2})
) where {N, N2, T <: Number}
    x_cur = x_el + u_el
    ndim = 1
    nnpe = length(x_cur) ÷ ndim
    h2 = typemax(T)
    for i in 1:(nnpe - 1)
        xi = x_cur[(i - 1) * ndim + 1]
        for j in (i + 1):nnpe
            dx = x_cur[(j - 1) * ndim + 1] - xi
            h2 = min(h2, dx * dx)
        end
    end
    return sqrt(h2)
end

function characteristic_element_length(
    x_el::SMatrix{2, N, T, N2},
    u_el::SMatrix{2, N, T, N2} = zero(u_el::SMatrix{2, N, T, N2})
) where {N, N2, T <: Number}
    x_cur = x_el + u_el
    ndim = 2
    nnpe = length(x_cur) ÷ ndim
    h2 = typemax(T)
    for i in 1:(nnpe - 1)
        xi = x_cur[(i - 1) * ndim + 1]
        yi = x_cur[(i - 1) * ndim + 2]
        for j in (i + 1):nnpe
            dx = x_cur[(j - 1) * ndim + 1] - xi
            dy = x_cur[(j - 1) * ndim + 2] - yi
            h2 = min(h2, dx * dx + dy * dy)
        end
    end
    return sqrt(h2)
end

function characteristic_element_length(
    x_el::SMatrix{3, N, T, N2},
    u_el::SMatrix{3, N, T, N2} = zero(u_el::SMatrix{3, N, T, N2})
) where {N, N2, T <: Number}
    x_cur = x_el + u_el
    ndim = 3
    nnpe = length(x_cur) ÷ ndim
    h2 = typemax(T)
    for i in 1:(nnpe - 1)
        xi = x_cur[(i - 1) * ndim + 1]
        yi = x_cur[(i - 1) * ndim + 2]
        zi = x_cur[(i - 1) * ndim + 3]
        for j in (i + 1):nnpe
            dx = x_cur[(j - 1) * ndim + 1] - xi
            dy = x_cur[(j - 1) * ndim + 2] - yi
            dz = x_cur[(j - 1) * ndim + 3] - zi
            h2 = min(h2, dx * dx + dy * dy + dz * dz)
        end
    end
    return sqrt(h2)
end

function setup_function(fspace::FunctionSpace, physics::AbstractPhysics)
    @assert false "setup_function needs to be implemented for $(typeof(physcs)). This method should return the functions necessary for this physics"
end



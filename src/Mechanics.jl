# function energy(cell::C, u_el) where C <: CellInterpolants

# end

# function energy(section::S, U::M) where {S <: TotalLagrangeSection, M <: AbstractMatrix}
#   # U_els = U[:, section.connectivity]
#   N, D = num_nodes(section), num_dimensions(section)
#   L    = N * D
#   for conn in section.connectivity
#     # unpack element level solution into static array
#     u_el = SMatrix{D, N, Float64, L}(@views U[:, conn])

#   end
# end
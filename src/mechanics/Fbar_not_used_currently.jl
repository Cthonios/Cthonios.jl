# function modify_deformation_gradient(::Type{<:FiniteElementContainers.AbstractMechanicsFormulation}, F_q, F_c)
#   if section.formulation <: FiniteElementContainers.PlaneStrain
#     F_q = ((det(F_c) / det(F_q))^(1. / 2.)) * F_q
#   else
#     @assert false "unsupported formulation"
#   end
#   return F_q
# end

function Fbar_deformation_gradient(::FiniteElementContainers.PlaneStrain, F_q, F_c)
  F_q = ((det(F_c) / det(F_q))^(1. / 2.)) * F_q
  return F_q
end

function Fbar_stifness(::FiniteElementContainers.PlaneStrain, G, G_c, P_v, A_m)
  T   = eltype(P_v)
  TxI = SMatrix{4, 4, T, 16}((
    P_v[1],  P_v[2],  P_v[3],  P_v[4],
    zero(T), zero(T), zero(T), zero(T),
    zero(T), zero(T), zero(T), zero(T),
    P_v[1],  P_v[2],  P_v[3],  P_v[4],
  ))
  T   = eltype(A_m)
  A_II = SMatrix{4, 4, T, 16}((
    A_m[1, 1] + A_m[1, 4], A_m[2, 1] + A_m[2, 4], A_m[3, 1] + A_m[3, 4], A_m[4, 1] + A_m[4, 4],
    zero(T),               zero(T),               zero(T),               zero(T),
    zero(T),               zero(T),               zero(T),               zero(T),
    A_m[1, 1] + A_m[1, 4], A_m[2, 1] + A_m[2, 4], A_m[3, 1] + A_m[3, 4], A_m[4, 1] + A_m[4, 4]
  ))
  Q = 0.5 * A_II - 0.5 * TxI

  K_bar = G * Q * (G_c - G)'
  return K_bar
end

function residual_and_stiffness(section::Section, X_el, U_el, q)
  
  # unpack some stuff, TODO state vars need to not be Hardcoded
  model, props = section.model, section.props
  state        = SVector{0, Float64}()

  # interpolants, also get centroid element
  ∇N_X, JxW, G   = interpolants(section.fspace.ref_fe, section.formulation, X_el, q)
  ∇N_X_c, _, G_c = interpolants(section.centroid_ref_fe, section.formulation, X_el, 1)

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
  ∇u_c = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X_c)
  F_q  = ∇u_q + one(∇u_q)
  F_c  = ∇u_c + one(∇u_c)

  # TODO wrap in if statement with Fbar as option
  # F_bar = F_q
  # F_bar = Fbar_deformation_gradient(section.formulation, F_q, F_c)
  F_bar = F_q

  # constitutive model
  P_q  = ConstitutiveModels.pk1_stress(model, props, F_bar, state)
  A_q  = Tensors.hessian(z -> ConstitutiveModels.energy(model, props, z, state), F_q)

  # vectorization
  P_v  = FiniteElementContainers.extract_stress(section.formulation, P_q)
  A_m  = FiniteElementContainers.extract_stiffness(section.formulation, A_q)

  # accumulation to element level
  R_q = JxW * G * P_v
  # K_q = JxW * G * A_m * G'

  K     = G * A_m * G'

  # TODO add if statement here for Fbar as an option
  K_bar = Fbar_stifness(section.formulation, G, G_c, P_v, A_m)

  # K_q   = JxW * (K + K_bar)
  K_q = JxW * K

  return R_q, K_q
end
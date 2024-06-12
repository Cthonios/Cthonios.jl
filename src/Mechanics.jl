# utils
function element_coordinates(section, X, conn)
  ND, NN, _, _ = size(section)
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[:, conn])
  return X_el
end

function element_fields(section, U, dof_conn)
  _, NN, _, _ = size(section)
  NF = num_fields(U)
  U_el = SMatrix{NF, NN, eltype(U), NF * NN}(@views U[dof_conn])
  return U_el
end

function element_properties(section, props, section_name, e)
  _, _, NP, _ = size(section)
  props_el = SVector{NP, eltype(props)}(@views props[section_name][:, e])
  return props_el
end

function interpolants(section, X_el, q)
  ref_fe, formulation = section.fspace.ref_fe, section.formulation
  N     = FiniteElementContainers.ReferenceFiniteElements.shape_function_values(ref_fe, q)
  ∇N_ξ  = FiniteElementContainers.ReferenceFiniteElements.shape_function_gradients(ref_fe, q)
  J     = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X  = (J_inv * ∇N_ξ')'
  JxW   = det(J) * FiniteElementContainers.ReferenceFiniteElements.quadrature_weights(ref_fe, q)
  G     = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return N, ∇N_X, JxW, G
end

# TODO maybe move somewhere else
# i.e. a solver? or linear solver?
function gradient end
function gradient! end
function internal_force end
function internal_force! end
function residual end
function residual! end

# Top level methods with domain exposed
function internal_energy!(domain::Domain, Uu::AbstractArray{Float64, 1})
  internal_energy!(domain.static, domain.cache, Uu)
  return nothing
end

function internal_force!(domain::Domain, Uu)
  internal_force!(domain.static, domain.cache, Uu)
  return nothing
end

function stiffness!(domain::Domain, Uu)
  stiffness!(domain.static, domain.cache, Uu)
  return nothing
end

function internal_energy(domain::Domain, Uu)
  internal_energy!(domain, Uu)
  return domain.cache.solver_cache.Π[1]
end

function residual(domain::Domain, Uu)
  internal_force!(domain, Uu)
  return @views domain.cache.solver_cache.assembler.residuals[domain.static.dof.unknown_dofs]
end

function stiffness(domain::Domain, Uu)
  stiffness!(domain, Uu)
  return SparseArrays.sparse!(domain.cache.solver_cache.assembler)
end

"""
useful for calculating residuals
"""
function internal_energy!(static, cache, Uu)
  update_fields!(static, cache, Uu)
  internal_energy!(static, cache)
  return nothing
end

function internal_force!(static, cache, Uu::Vector{Float64})
  update_fields!(static, cache, Uu)
  internal_force!(static, cache)
  return nothing
end

function stiffness!(static, cache, Uu)
  update_fields!(static, cache, Uu)
  stiffness!(static, cache)
  return nothing
end

# function stiffness_action!(domain::Domain, Uu)

# end

# section iterators
function internal_energy!(static, cache)
  cache.solver_cache.Π .= zero(eltype(cache.solver_cache.Π))
  for (section_name, section) in pairs(static.sections)
    _, _, NP, NS = size(section)
    fspace = section.fspace

    for e in 1:num_elements(fspace)
      conn = connectivity(fspace, e)
      dof_conn = dof_connectivity(fspace, e)
      X = element_coordinates(section, cache.X, conn)
      U = element_fields(section, cache.U, dof_conn)
      # @time props_el = element_properties(section, cache.props, section_name, e)
      # below allocates if we wrap it in a method?
      props_el = SVector{NP, eltype(cache.props)}(@views cache.props[section_name][:, e])
      for q in 1:num_q_points(fspace)
        state_old_q = SVector{NS, eltype(cache.state_old)}(@views cache.state_old[section_name][:, q, e])
        # get inerpolants
        N, ∇N_X, JxW, G = interpolants(section, X, q)
        # calculate values at quadrature point
        X_q  = X * N
        u_q  = U * N
        ∇u_q = U * ∇N_X
        # run quadrature level routine
        W_q, state_new_q = strain_energy(
          section,
          cache.times.Δt, X_q, u_q, ∇u_q, props_el, state_old_q
        )
        @views cache.solver_cache.Πs[section_name][q, e] = JxW * W_q
        @views cache.state_new[section_name][:, q, e] .= state_new_q
      end
    end
  end
  cache.solver_cache.Π[1] = sum(cache.solver_cache.Πs)
  return nothing
end

function internal_force!(static, cache)
  cache.solver_cache.assembler.residuals .= 
    zero(eltype(cache.solver_cache.assembler.residuals))
  for (section_name, section) in pairs(static.sections)
    ND, NN, NP, NS = size(section)
    NF = num_fields(cache.U)
    fspace = section.fspace

    for e in 1:num_elements(fspace)
      conn = connectivity(fspace, e)
      dof_conn = dof_connectivity(fspace, e)
      X = element_coordinates(section, cache.X, conn)
      U = element_fields(section, cache.U, dof_conn)
      # below allocates if we wrap it in a method?
      props_el = SVector{NP, eltype(cache.props)}(@views cache.props[section_name][:, e])
      f_el = zeros(SVector{NF * NN, eltype(cache.solver_cache.assembler.residuals)})
      for q in 1:num_q_points(fspace)
        state_old_q = SVector{NS, eltype(cache.state_old)}(@views cache.state_old[section_name][:, q, e])
        
        # get inerpolants
        N, ∇N_X, JxW, G = interpolants(section, X, q)
        # calculate values at quadrature point
        X_q  = X * N
        u_q  = U * N
        ∇u_q = U * ∇N_X
        # run quadrature level routine
        P_q, state_new_q = pk1_stress(
          section,
          cache.times.Δt, X_q, u_q, ∇u_q, props_el, state_old_q
        )
        P_v = FiniteElementContainers.extract_stress(section.formulation, P_q)
        f_q = G * P_v
        f_el = f_el + JxW * f_q

        @views cache.state_new[section_name][:, q, e] .= state_new_q
      end
      assemble!(cache.solver_cache.assembler.residuals, f_el, dof_conn)
    end
  end
  return nothing
end

function stiffness!(static, cache)
  cache.solver_cache.assembler.stiffnesses .= 
    zero(eltype(cache.solver_cache.assembler.stiffnesses))
  for (block_num, (section_name, section)) in enumerate(pairs(static.sections))
    ND, NN, NP, NS = size(section)
    NF = num_fields(cache.U)
    fspace = section.fspace

    for e in 1:num_elements(fspace)
      conn = connectivity(fspace, e)
      dof_conn = dof_connectivity(fspace, e)
      X = element_coordinates(section, cache.X, conn)
      U = element_fields(section, cache.U, dof_conn)
      # below allocates if we wrap it in a method?
      props_el = SVector{NP, eltype(cache.props)}(@views cache.props[section_name][:, e])
      K_el = zeros(SMatrix{NF * NN, NF * NN, eltype(cache.solver_cache.assembler.residuals), NF * NN * NF * NN})
      for q in 1:num_q_points(fspace)
        state_old_q = SVector{NS, eltype(cache.state_old)}(@views cache.state_old[section_name][:, q, e])
        
        # get inerpolants
        N, ∇N_X, JxW, G = interpolants(section, X, q)
        # calculate values at quadrature point
        X_q  = X * N
        u_q  = U * N
        ∇u_q = U * ∇N_X
        # run quadrature level routine
        K_q, state_new_q = material_tangent(
          section,
          cache.times.Δt, X_q, u_q, ∇u_q, props_el, state_old_q
        )
        K_v = FiniteElementContainers.extract_stiffness(section.formulation, K_q)
        K_q = G * K_v * G'
        K_el = K_el + JxW * K_q

        @views cache.state_new[section_name][:, q, e] .= state_new_q
      end
      assemble!(cache.solver_cache.assembler, K_el, block_num, e)
    end
  end
  return nothing
end

# quadrature point kernels
function strain_energy(
  section::TotalLagrangeSectionInternal,
  Δt, X_q, u_q, ∇u_q, props_el, state_old_q
)
  # TODO temp hardcoded temperature
  # Hook up stream with everything elase
  θ_q = 0.0
  # kinematics
  ∇u_q = modify_field_gradients(section, ∇u_q)
  F_q = ∇u_q + one(∇u_q)
  # constitutive
  ψ_q, state_new_q = ConstitutiveModels.helmholtz_free_energy(
    section.model, props_el, Δt, F_q, θ_q, state_old_q
  )

  return ψ_q, state_new_q
end

function pk1_stress(
  section::TotalLagrangeSectionInternal,
  Δt, X_q, u_q, ∇u_q, props_el, state_old_q
)
  # TODO temp hardcoded time step and temperature
  # Hook up stream with everything elase
  θ_q = 0.0
  # kinematics
  ∇u_q = modify_field_gradients(section, ∇u_q)
  F_q = ∇u_q + one(∇u_q)
  # constitutive
  P_q, state_new_q = ConstitutiveModels.pk1_stress(
    section.model, props_el, Δt, F_q, θ_q, state_old_q
  )
  return P_q, state_new_q
end

function material_tangent(
  section::TotalLagrangeSectionInternal,
  Δt, X_q, u_q, ∇u_q, props_el, state_old_q
)
  # TODO temp hardcoded time step and temperature
  # Hook up stream with everything elase
  θ_q = 0.0
  # kinematics
  ∇u_q = modify_field_gradients(section, ∇u_q)
  F_q = ∇u_q + one(∇u_q)
  # constitutive
  K_q, state_new_q = ConstitutiveModels.material_tangent(
    section.model, props_el, Δt, F_q, θ_q, state_old_q
  )
  return K_q, state_new_q
end

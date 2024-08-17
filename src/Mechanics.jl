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

# defined in extension
function gradient end
function gradient! end
function hvp! end
function hvp end
function internal_force end
function residual! end

function grad end

# Out of place methods with domain and Uu exposed
for (op1, op2, op3) in zip(
  (:internal_energy, :residual, :stiffness),
  (:internal_energy!, :internal_force!, :stiffness!),
  (
    Expr(:block, [
      Meta.parse("domain.solver_cache.Π[1] = sum(domain.solver_cache.Πs)"),
      Meta.parse("return domain.solver_cache.Π[1]")
    ]...),
    Meta.parse("@views domain.solver_cache.assembler.residuals[domain.dof.unknown_dofs]"),
    Meta.parse("SparseArrays.sparse!(domain.solver_cache)")
  )
) 
  @eval begin
    function ($op1)(domain, Uu)
      ($op2)(domain, Uu)
      ($op3)
    end
  end
end

# In place methods with domain and Uu exposed
for (op1, op2) in zip(
  (:internal_energy!, :internal_force!, :stiffness!),
  (
    Meta.parse("domain.solver_cache.Π[1] = sum(domain.solver_cache.Πs)"),
    Meta.parse(""), Meta.parse("")
  )
)
  @eval begin 
    function ($op1)(domain::Domain, Uu) 
      update_fields!(domain, Uu)
      ($op1)(domain)
      ($op2)
    end
  end
end

# note V is already a field here
function stiffness_action!(domain::Domain, Vv)
  # update_fields!(domain, Uu)
  FiniteElementContainers.update_fields!(domain.solver_cache.V, domain.dof, Vv)
  stiffness_action!(domain)
  return nothing
end

# section iterators
for (op1, op2, op3, op4, op5, op6) in zip(
  ( # name of method to define
    :internal_energy!,
    :internal_force!,
    :stiffness!
  ),
  ( # reset operation
    Meta.parse("domain.solver_cache.Π .= zero(eltype(domain.solver_cache.Π))"),
    Meta.parse("domain.solver_cache.assembler.residuals .= zero(eltype(domain.solver_cache.assembler.residuals))"),
    Meta.parse("domain.solver_cache.assembler.stiffnesses .= zero(eltype(domain.solver_cache.assembler.stiffnesses))")
  ),
  (
    # optional scratch var
    Meta.parse(""),
    Meta.parse("f_el = zeros(SVector{NF * NN, eltype(domain.solver_cache.assembler.residuals)})"),
    Meta.parse("K_el = zeros(SMatrix{NF * NN, NF * NN, eltype(domain.solver_cache.assembler.stiffnesses), NF * NN * NF * NN})")
  ),
  ( # name of quadrature point method to use
    :strain_energy,
    :pk1_stress,
    :material_tangent
  ),
  (
    # optional quadrature point clean up
    Meta.parse("@views domain.solver_cache.Πs[section_name][q, e] = JxW * qoi_q"),
    Expr(:block, [
      Meta.parse("P_v = FiniteElementContainers.extract_stress(section.formulation, qoi_q)"),
      Meta.parse("P_q = G * P_v"),
      Meta.parse("f_el = f_el + JxW * P_q")
    ]...),
    Expr(:block, [
      Meta.parse("K_v = FiniteElementContainers.extract_stiffness(section.formulation, qoi_q)")
      Meta.parse("K_q = G * K_v * G'")
      Meta.parse("K_el = K_el + JxW * K_q")
    ]...)
  ),
  (
    # optional element clean up
    Meta.parse(""),
    Meta.parse("assemble!(domain.solver_cache.assembler.residuals, f_el, dof_conn)"),
    Meta.parse("assemble!(domain.solver_cache.assembler, K_el, block_num, e)")
  ),
)
  @eval begin
    function ($op1)(domain)
      # reset
      ($op2)

      # loop over sections
      for (block_num, (section_name, section)) in enumerate(pairs(domain.sections))
        ND, NN, NP, NS = size(section)
        NF = num_fields(domain.U)
        fspace = section.fspace
    
        for e in 1:num_elements(fspace)
          conn = connectivity(fspace, e)
          dof_conn = dof_connectivity(fspace, e)
          X = element_coordinates(section, domain.X, conn)
          U = element_fields(section, domain.U, dof_conn)
          # below allocates if we wrap it in a method?
          props_el = SVector{NP, eltype(domain.props)}(@views domain.props[section_name][:, e])
          # props_el = SVector{2, Float64}((0.833, 0.3846))
          # temporary storage
          ($op3)
          for q in 1:num_q_points(fspace)
            # state_old_q = SVector{NS, eltype(domain.state_old)}(@views domain.state_old[section_name][:, q, e])
            state_old_q = SVector{0, Float64}()
            # get inerpolants
            N, ∇N_X, JxW, G = interpolants(section, X, q)
            # calculate values at quadrature point
            X_q  = X * N
            u_q  = U * N
            ∇u_q = U * ∇N_X
            # run quadrature level routine
            qoi_q, state_new_q = ($op4)(
              section,
              domain.times.Δt, X_q, u_q, ∇u_q, props_el, state_old_q
            )

            # set values
            ($op5)

            # common across all methods
            # @views domain.state_new[section_name][:, q, e] .= state_new_q
          end
          ($op6)
        end
      end
      return nothing
    end
  end
end

function stiffness_action!(domain::Domain)
  domain.solver_cache.Hv .= zero(eltype(domain.solver_cache.Hv))

  # loop over sections
  for (block_num, (section_name, section)) in enumerate(pairs(domain.sections))
    ND, NN, NP, NS = size(section)
    NF = num_fields(domain.U)
    fspace = section.fspace

    for e in 1:num_elements(fspace)
      conn = connectivity(fspace, e)
      dof_conn = dof_connectivity(fspace, e)
      X = element_coordinates(section, domain.X, conn)
      U = element_fields(section, domain.U, dof_conn)
      # V = element_fields(section, domain.solver_cache.V, dof_conn)
      V = SVector{NF * NN, eltype(domain.solver_cache.V)}(@views domain.solver_cache.V[dof_conn])
      # below allocates if we wrap it in a method?
      props_el = SVector{NP, eltype(domain.props)}(@views domain.props[section_name][:, e])
      # props_el = SVector{2, Float64}((0.833, 0.3846))
      # temporary storage
      # ($op3)
      Hv_el = SVector{NF * NN, eltype(domain.solver_cache.Hv)}(@views domain.solver_cache.Hv[dof_conn])
      for q in 1:num_q_points(fspace)
        # state_old_q = SVector{NS, eltype(domain.state_old)}(@views domain.state_old[section_name][:, q, e])
        state_old_q = SVector{0, Float64}()
        # get inerpolants
        N, ∇N_X, JxW, G = interpolants(section, X, q)
        # calculate values at quadrature point
        X_q  = X * N
        u_q  = U * N
        ∇u_q = U * ∇N_X
        # run quadrature level routine
        A_q, state_new_q = material_tangent(
          section,
          domain.times.Δt, X_q, u_q, ∇u_q, props_el, state_old_q
        )
        A_v = FiniteElementContainers.extract_stiffness(section.formulation, A_q)
        A_q = G * A_v * G' * V
        Hv_el = Hv_el + JxW * A_q

        # common across all methods
        # @views domain.state_new[section_name][:, q, e] .= state_new_q
      end
      assemble!(domain.solver_cache.Hv, Hv_el, dof_conn)
    end
  end
  return nothing
end

# quadrature point kernels
for (op1, op2) in zip(
  (
    :strain_energy, 
    :pk1_stress, 
    :material_tangent
  ),
  (
    :(ConstitutiveModels.helmholtz_free_energy),
    :(ConstitutiveModels.pk1_stress),
    :(ConstitutiveModels.material_tangent)
  )
)
  @eval begin
    function ($op1)(
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
      qoi_q, state_new_q = ($op2)(
        section.model, props_el, Δt, F_q, θ_q, state_old_q
      )
      return qoi_q, state_new_q
    end
  end
end

# top level in place method that loops over sections
# function energy!(Π, U, domain::QuasiStaticDomain, Uu, state, props, X)
function energy_old!(Πs, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  Πs .= zero(eltype(Πs))
  for (name, section) in pairs(sections)
    Πs_temp        = Πs[name]
    state_old_temp = state_old[name]
    state_new_temp = state_new[name]
    props_temp     = props[name]
    energy!(
      Πs_temp, state_new_temp, 
      section, Δt, X, U, props_temp, state_old_temp, backend
    )

    state_new[name] .= state_new_temp
    Πs[name] .= Πs_temp
  end 
  return nothing
end 

function energy!(Πs, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  Πs .= zero(eltype(Πs))
  for (name, section) in pairs(sections)
    ND, NN, NP, NS = size(section)

    for e in 1:num_elements(section)
      conn = dof_connectivity(section, e)
      X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
      U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
      props_el = SVector{NP, eltype(props)}(@views props[name][:, e])
      for q in 1:FiniteElementContainers.num_q_points(section)
        state_old_q = SVector{NS, eltype(state_old)}(@views state_old[name][:, q, e])
        W_q, state_new_q = energy(
          section.model, section.formulation,
          Δt, X_el, U_el, props_el, state_old_q,
          interpolants(section, q)...
        )
        @views Πs[name][q, e] = W_q
        @views state_new[name][:, q, e] .= state_new_q
      end
    end
  end 
  return nothing
end 

function internal_force!(f, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  f .= zero(eltype(f))
  for (name, section) in pairs(sections)
    ND, NN, NP, NS = size(section)

    for e in 1:num_elements(section)
      conn = dof_connectivity(section, e)
      X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
      U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
      props_el = SVector{NP, eltype(props)}(@views props[name][:, e])

      f_el = zeros(SVector{ND * NN, eltype(f)})

      for q in 1:FiniteElementContainers.num_q_points(section)
        state_old_q = SVector{NS, eltype(state_old)}(@views state_old[name][:, q, e])
        f_q, state_new_q = internal_force(
          section.model, section.formulation,
          Δt, X_el, U_el, props_el, state_old_q,
          interpolants(section, q)...
        )
        f_el = f_el + f_q
        @views state_new[name][:, q, e] .= state_new_q
      end

      # assemble into global field
      assemble!(f, f_el, conn)
    end
  end 
  return nothing
end 

function stiffness!(assembler::StaticAssembler, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)

    ND, NN, NP, NS = size(section)
    NDOF = ND * NN

    for e in 1:num_elements(section)
      conn = dof_connectivity(section, e)
      X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
      U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
      props_el = SVector{NP, eltype(props)}(@views props[name][:, e])

      K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

      for q in 1:FiniteElementContainers.num_q_points(section)
        state_old_q = SVector{NS, eltype(state_old)}(@views state_old[name][:, q, e])

        K_q, state_new_q = stiffness(
          section.model, section.formulation,
          Δt, X_el, U_el, props_el, state_old_q,
          interpolants(section, q)...
        )
        K_el = K_el + K_q
        @views state_new[name][:, q, e] .= state_new_q
      end

      # assemble into global matrix
      assemble!(assembler, K_el, block_count, e)
    end 

    block_count = block_count + 1
  end 
  return nothing
end 

function stiffness_action!(Kv, state_new, sections, Δt, X, U, props, state_old, V, backend::Backend)
  Kv .= zero(eltype(Kv))
  for (name, section) in pairs(sections)

    ND, NN, NP, NS = size(section)
    NDOF = ND * NN

    for e in 1:num_elements(section)
      conn = dof_connectivity(section, e)
      X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
      U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
      V_el = SVector{ND * NN, eltype(V)}(@views V[conn])
      props_el = SVector{NP, eltype(props)}(@views props[name][:, e])

      Kv_el = zeros(SVector{NDOF, eltype(U)})

      for q in 1:FiniteElementContainers.num_q_points(section)
        state_old_q = SVector{NS, eltype(state_old)}(@views state_old[name][:, q, e])
        
        K_q, state_new_q = stiffness(
          section.model, section.formulation,
          Δt, X_el, U_el, props_el, state_old_q,
          interpolants(section, q)...
        )
        Kv_el = Kv_el + K_q * V_el
        @views state_new[name][:, q, e] .= state_new_q
      end

      # assemble into global matrix
      assemble!(Kv, Kv_el, conn)
    end
  end 
  return nothing
end

# dual return outputs
function energy_and_internal_force!(Πs, f, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  Πs .= zero(eltype(Πs))
  f .= zero(eltype(f))
  for (name, section) in pairs(sections)

    ND, NN, NP, NS = size(section)
    NDOF = ND * NN

    for e in 1:num_elements(section)
      conn = dof_connectivity(section, e)
      X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
      U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
      props_el = SVector{NP, eltype(props)}(@views props[name][:, e])

      f_el = zeros(SVector{NDOF, eltype(U)})

      for q in 1:FiniteElementContainers.num_q_points(section)
        state_old_q = SVector{NS, eltype(state_old)}(@views state_old[name][:, q, e])
        W_q, f_q, state_new_q = energy_and_internal_force(
          section.model, section.formulation,
          Δt, X_el, U_el, props_el, state_old_q,
          interpolants(section, q)...
        )
        @views Πs[name][q, e] = W_q
        @views state_new[name][:, q, e] .= state_new_q
        f_el = f_el + f_q
      end

      # assemble into global matrix
      assemble!(f, f_el, conn)
    end
  end
  return nothing
end

function internal_force_and_stiffness!(f, assembler, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  f .= zero(eltype(f))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    ND, NN, NP, NS = size(section)
    NDOF = ND * NN

    ND, NN, NP, NS = size(section)

    for e in 1:num_elements(section)
      conn = dof_connectivity(section, e)
      X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
      U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
      props_el = SVector{NP, eltype(props)}(@views props[name][:, e])

      f_el = zeros(SVector{NDOF, eltype(f)})
      K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

      for q in 1:FiniteElementContainers.num_q_points(section)
        state_old_q = SVector{NS, eltype(state_old)}(@views state_old[name][:, q, e])
        f_q, K_q, state_new_q = internal_force_and_stiffness(
          section.model, section.formulation,
          Δt, X_el, U_el, props_el, state_old_q,
          interpolants(section, q)...
        )
        f_el = f_el + f_q
        K_el = K_el + K_q
        @views state_new[name][:, q, e] .= state_new_q
      end

      # assemble into global fields
      assemble!(f, f_el, conn)
      assemble!(assembler, K_el, block_count, e)
    end

    block_count = block_count + 1
  end 
  return nothing
end 

# triple return outputs
function energy_internal_force_and_stiffness!(Πs, f, assembler, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  Πs .= zero(eltype(Πs))
  f .= zero(eltype(f))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    ND, NN, NP, NS = size(section)
    NDOF = ND * NN

    for e in 1:num_elements(section)
      conn = dof_connectivity(section, e)
      X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
      U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
      props_el = SVector{NP, eltype(props)}(@views props[name][:, e])

      f_el = zeros(SVector{NDOF, eltype(f)})
      K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})
  
      for q in 1:FiniteElementContainers.num_q_points(section)
        state_old_q = SVector{NS, eltype(state_old)}(@views state_old[name][:, q, e])
        W_q, f_q, K_q, state_new_q = energy_internal_force_and_stiffness(
          section.model, section.formulation,
          Δt, X_el, U_el, props_el, state_old_q,
          interpolants(section, q)...
        )
        @views Πs[name][q, e] = W_q
        f_el = f_el + f_q
        K_el = K_el + K_q
        @views state_new[name][:, q, e] .= state_new_q
      end

      # assemble into global fields
      assemble!(f, f_el, conn)
      assemble!(assembler, K_el, block_count, e)
    end

    block_count = block_count + 1
  end
  return nothing
end

function energy_internal_force_and_stiffness_action!(Πs, f, Hv, state_new, sections, Δt, X, U, props, state_old, V, backend::Backend)
  Πs .= zero(eltype(Πs))
  f .= zero(eltype(f))
  Hv .= zero(eltype(Hv))
  for (name, section) in pairs(sections)
    ND, NN, NP, NS = size(section)
    NDOF = ND * NN

    for e in 1:num_elements(section)
      conn = dof_connectivity(section, e)
      X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
      U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
      V_el = SVector{NDOF, eltype(V)}(@views V[conn])
      props_el = SVector{NP, eltype(props)}(@views props[name][:, e])

      f_el = zeros(SVector{NDOF, eltype(U)})
      Hv_el = zeros(SVector{NDOF, eltype(U)})

      for q in 1:FiniteElementContainers.num_q_points(section)
        state_old_q = SVector{NS, eltype(state_old)}(@views state_old[name][:, q, e])
        W_q, f_q, H_q, state_new_q = energy_internal_force_and_stiffness(
          section.model, section.formulation,
          Δt, X_el, U_el, props_el, state_old_q,
          interpolants(section, q)...
        )
        Πs[name][q, e] = W_q
        f_el = f_el + f_q
        Hv_el = Hv_el + H_q * V_el
        state_new[name][:, q, e] .= state_new_q
      end
    end
  end
  return nothing
end

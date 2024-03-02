@kernel function energy_kernel!(
  Πs, 
  section::TotalLagrangeSection,
  U, state, props, X,
)
  # indices
  q, e = @index(Global, NTuple)

  # unpack arrays
  _, U_el, state_q, props_el, X_el = unpack(section, U, state, props, X, q, e)

  # get interpolants
  N, ∇N_ξ, w = interpolants(section, q)

  # run routine
  Π_q = energy(
    section.model, section.formulation,
    U_el, state_q, props_el, X_el,
    N, ∇N_ξ, w
  )

  # update arrays
  Πs[q, e] = Π_q

  # return nothing # issue in kernel abstractions
end

function energy!(Π, Πs, sections, U, state, props, X, backend::CPU)

  # kernel setup
  kernel! = energy_kernel!(backend)

  for (name, section) in pairs(sections)
    NQ = FiniteElementContainers.num_q_points(section)
    NE = num_elements(section)
    Πs_temp    = @views Πs[name]
    state_temp = @views state[name]
    props_temp = @views props[name]
    kernel!(
      Πs_temp, 
      section,
      U, state_temp, props_temp, X,
      ndrange=(NQ, NE)
    )
  end 
  synchronize(backend)
  Π[1] = sum(Πs)
  return nothing
end 

@kernel function energy_and_internal_force_kernel!(
  Πs, f,
  section::TotalLagrangeSection,
  U, state, props, X,
)
  # indices
  q, e = @index(Global, NTuple)

  # unpack arrays
  conn, U_el, state_q, props_el, X_el = unpack(section, U, state, props, X, q, e)
  # conn = conn |> Vector
  # get interpolants
  N, ∇N_ξ, w = interpolants(section, q)

  # run routine
  Π_q, f_q = energy_and_internal_force(
    section.model, section.formulation,
    U_el, state_q, props_el, X_el,
    N, ∇N_ξ, w
  )

  conn = conn |> Vector
  f_q = f_q |> Vector
  # @show typeof(f_q)
  # @show typeof(f)
  # vals = f.vals
  # update arrays
  Πs[q, e] = Π_q
  
  for n in 1:size(conn, 1)
    # m = conn[n]
    Atomix.@atomic f[conn[n]] += f_q[n]
  end

  # Atomix.@atomic f.vals[conn] += f_q
  # return nothing # issue in kernel abstractions
end

function energy_and_internal_force!(Πs, f, sections, U, state, props, X, backend::CPU)

  # kernel setup
  kernel! = energy_and_internal_force_kernel!(backend)

  for (name, section) in pairs(sections)
    NQ = FiniteElementContainers.num_q_points(section)
    NE = num_elements(section)
    Πs_temp    = @views Πs[name]
    state_temp = @views state[name]
    props_temp = @views props[name]
    kernel!(
      Πs_temp, f,
      section,
      U, state_temp, props_temp, X,
      ndrange=(NQ, NE)
    )
  end 
  synchronize(backend)
  # Π[1] = sum(Πs)
  return nothing
end 
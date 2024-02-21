# for kernel abstractions WIP
function energy(domain::QuasiStaticDomain, backend)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  energy(domain, Uu, state, props, X, backend)
end   

# for kernel abstractions WIP
function energy(domain::QuasiStaticDomain, Uu::V, state, props, X, backend) where V <: AbstractVector
  U = create_fields(domain)
  Π = zeros(eltype(Uu), 1)
  update_fields!(U, domain, X, Uu)

  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy!(Π, section, U, state_temp, props_temp, X, backend)
  end 
end

# for kernel abstractions WIP
function energy!(Π, U, domain::QuasiStaticDomain, Uu, state, props, X, backend)
  update_fields!(U, domain, X, Uu)
  Π .= zero(eltype(Π))
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy!(Π, section, U, state_temp, props_temp, X, backend)
  end 
  return nothing
end 


# for kernel abstractions WIP
@kernel function energy_kernel!(energies, section, U, state, props, X)
  Q, E = @index(Global, NTuple)

  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)

  conn = dof_connectivity(section, E)

  U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
  props_el = SVector{NP, eltype(props)}(@views props[:, E])
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
  state_q = SVector{NS, eltype(props)}(@views state[:, Q, E])

  Atomix.@atomic energies[1] = energy(
    section.model, section.formulation,
    U_el, state_q, props_el, X_el,
    shape_function_values(section, Q),
    shape_function_gradients(section, Q),
    quadrature_weights(section, Q)
  )

  return nothing
end

# for kernel abstractions WIP
function energy!(Π, section::TotalLagrangeSection, U, state, props, X, backend)
  kernel! = energy_kernel!(backend)
  kernel!(Π, section, U, state, props, X, ndrange=(FiniteElementContainers.num_q_points(section), num_elements(section)))
  return nothing
end

# other GPU attempts
function energy_mapreduce(section::TotalLagrangeSection, U, state, props, X)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)

  conn = dof_connectivity(section)

  U_els = reinterpret(SMatrix{ND, NN, eltype(U), ND * NN}, @views vec(U[conn]))
  X_els = reinterpret(SMatrix{ND, NN, eltype(U), ND * NN}, @views vec(X[conn]))
  # state = reinterpret(SVector{})

  model = section.model
  form = section.formulation

  Ns = shape_function_values(section)
  ∇N_ξs = shape_function_gradients(section)
  ws = quadrature_weights(section)

  # func = (a, b, c, d, e, f, g, h, i) -> map(energy, (a,), (b,), (c,), eachcol(d), (e,), (f,), g, h, i)
  # W = map(
  #   func,
  #   # (a, b, c, d, e, f, g, h, i) -> map(energy, (a,), (b,), (c,), d, (e,), (f,), g, h, i),
  #   (model,), (form,),
  #   U_els, eachslice(state; dims=(3,)), eachcol(props), X_els,
  #   (Ns,), (∇N_ξs,), (ws,)
  # )

  # func = (a, )
  # W = map((x, y) -> (x, y), (model,), U_els)
  
  # W = mapreduce(
  #   (u, s, p, x) -> mapreduce(
  #     (s_q, N, ∇N, w) -> energy(model, form, u, s_q, p, x, N, ∇N, w), 
  #     +,
  #     eachcol(s), Ns, ∇N_ξs, ws
  #   ), +,
  #   U_els, eachslice(state; dims=(3,)), eachcol(props), X_els
  # )
  

  # s_temp = state[:, :, 1]
  # p = props[:, 1]
  # u = U_els[1]
  # x_temp = X_els[1]

  # # inner_func = (f, op, s, Ns, ∇Ns, ws) -> mapreduce(f, op, eachcol(s), Ns, ∇Ns, ws)
  # # TODO need to patch up below
  # inner_func = (f, op, x) -> mapreduce(f, op, x)

  # # @show s_temp
  # @time W = inner_func(
  #   z -> energy(model, form, u, SVector{0, Float64}(), p, x, z[1], z[2], z[3]), +,
  #   zip(Ns, ∇N_ξs, ws)
  # )

  # W = mapreduce(
  #   y -> energy(model, form, u, SVector{0, Float64}, p, x, y[1], y[2], y[3]), +,
  #   zip(Ns, ∇N_ξs, ws)
  # )

  W = mapreduce(
    x -> mapreduce(
      y -> energy(model, form, x[1], SVector{0, Float64}, x[2], x[3], y[1], y[2], y[3]), +,
      zip(Ns, ∇N_ξs, ws)
    ), +,
    zip(U_els, eachcol(props), X_els)
  )

  # W = 0.0

  return W
end

function inner_mapreduce(f, op, s, Ns, ∇Ns, ws)
  return mapreduce(f, op, eachcol(s), Ns, ∇Ns, ws)
end

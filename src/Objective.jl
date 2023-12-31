# FunctionWrappers probably won't work except for maybe with enzyme?
struct Objective{Domain, Backend, F <: Function}
  domain::Domain
  backend::Backend
  func::F
end

objective(o::Objective, u) = o.func(o.domain, u)
objective(o::Objective, x, u) = o.func(o.domain, x, u)

grad_x(o::Objective, x, u) = AD.gradient(o.backend, z -> o.func(o.domain, z, u), x)[1]
grad_u(o::Objective, x, u) = AD.gradient(o.backend, z -> o.func(o.domain, x, z), u)[1]

# hvp_x is problem currentl
hvp_x(o::Objective, x, u) = AD.pushforward_function(o.backend, z -> grad_x(o, z, u), x)
hvp_u(o::Objective, x, u) = AD.pushforward_function(o.backend, z -> grad_u(o, x, z), u)

# func_map(o::Objective, x, u, m) = LinearMaps.FunctionMap(hvp_u(o, x, u), m)
# func_map(o::Objective, x, u) = LinearProblem(hvp_u(o, x, u), u)

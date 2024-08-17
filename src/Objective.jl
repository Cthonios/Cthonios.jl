# abstract type AbstractObjective{F1, F2, F3} end

struct Objective{F1, F2, F3} #<: AbstractObjective{F1, F2, F3}
  value::F1
  gradient::F2
  hessian::F3
end

# a little dirty but whatever
gradient!(obj::Objective, domain::Domain) = obj.gradient(domain)
hessian!(obj::Objective, domain::Domain) = obj.hessian(domain)
value!(obj::Objective, domain::Domain) = obj.value(domain)
gradient!(obj::Objective, domain::Domain, Uu) = obj.gradient(domain, Uu)
hessian!(obj::Objective, domain::Domain, Uu) = obj.hessian(domain, Uu)
value!(obj::Objective, domain::Domain, Uu) = obj.value(domain, Uu)

gradient(obj::Objective, domain::Domain, Uu) = obj.gradient(domain, Uu)
hessian(obj::Objective, domain::Domain, Uu) = obj.hessian(domain, Uu)
value(obj::Objective, domain::Domain, Uu) = obj.value(domain, Uu)

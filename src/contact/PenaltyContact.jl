struct PenaltyContact{E}
  edges::E
end

function PenaltyContact(d::Domain)
  return PenaltyContact(edges(d))
end

function constraint(::PenaltyContact, level_set, domain, Uu, p)
  # for each edge
  # evaluate level set on each quadrature point on the edge
  # so we'll need to have an edge ref_fe that maps reference
  # q points to CURRENT coordinates
  
end

struct WarmStart{U1, U2}
  Uu_old::U1
  U_old::U2
end

function WarmStart(o::Objective)
  return WarmStart(create_unknowns(o.domain), create_fields(o.domain))
end

# solve! in EnzymeExt currently.
# TODO figure out a way to do it analytically

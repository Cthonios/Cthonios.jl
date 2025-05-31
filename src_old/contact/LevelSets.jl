abstract type AbstractLevelSet end

struct Corner{T} <: AbstractLevelSet
end

function (::Corner{T})(x, args...) where T
  x_loc, y_loc = args
  return min(x[1] - x_loc, x[2] - y_loc)
end

struct Plane{T} <: AbstractLevelSet
  height::T
end

function (plane::Plane{T})(x, args...) where T
  return plane.height - x[2]
end

struct Sphere{T} <: AbstractLevelSet
  radius::T
end

function (sphere::Sphere{T})(x, args...) where T
  R = sphere.radius
  x_loc, y_loc = args
  r = sqrt((x[1] - x_loc)^2 + (x[2] - y_loc)^2) - R
  return r
end

struct Combined{L1, L2} <: AbstractLevelSet
  l1::L1
  l2::L2
end

function (combined::Combined{L1, L2})(x, args...) where {L1, L2}
  return min(combined.l1(x, args...), combined.l2(x, args...))
end

# containers for running
struct LevelSetContactPair{S, L <: AbstractLevelSet} <: AbstractContactInput
  sset_name::S
  l_set::L
end

struct LevelSetContactPairInternal{S, L <: AbstractLevelSet} <: AbstractContactInternal
  c_sec::S
  l_set::L
end

function LevelSetContactPairInternal(mesh, dof, section, contact_pair)
  sset = ContactSet(mesh, String(contact_pair.sset_name))
  sec_sset = SurfaceSectionInternal(mesh, dof, section, sset)
  return LevelSetContactPairInternal(sec_sset, contact_pair.l_set)
end

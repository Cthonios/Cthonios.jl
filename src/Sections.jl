# SECTION FORMULATIONS

abstract type AbstractSectionFormulation end

struct DefaultFormulation <: AbstractSectionFormulation
end

modify_field_gradient(::DefaultFormulation, ∇u_q) = ∇u_q

struct PlaneStrain <: AbstractSectionFormulation
end

# modify_field_gradient(::PlaneStrain, ∇u_q) = PaddedView(0.0, ∇u_q, (3, 3))

function modify_field_gradient(::PlaneStrain, ∇u_q::Tensor{2, 2, T, 4}) where T <: Number
  return Tensor{2, 3, T, 9}((
    ∇u_q[1, 1], ∇u_q[2, 1], 0.0,
    ∇u_q[1, 2], ∇u_q[2, 2], 0.0,
    0.0,        0.0,        0.0
  ))
end

###############################################################################

# Section values
abstract type AbstractSectionValues{NState, T} end
abstract type AbstractSectionQuadraturePointValues{NState, T} end

# TODO add state variables
struct SectionQuadraturePointValues{NState, T <: Number} <: AbstractSectionQuadraturePointValues{NState, T}
  W::T
  P::Tensor{2, 3, T, 9}
  A::Tensor{4, 3, T, 81}
  state_old::SVector{NState, T}
  state_new::SVector{NState, T}
end

struct PreAllocatedTotalLagrangeSectionValues{NState, T <: Number, S} <: AbstractSectionValues{NState, T}
  values::S
end

function Base.show(io::IO, ::PreAllocatedTotalLagrangeSectionValues)
  println(io, "PreAllocatedTotalLagrangeSectionValues")
end

# TODO maybe add different options to allow for different state 
# variable states, i.e. initial conditions
function PreAllocatedTotalLagrangeSectionValues(
  material_state::V, Q::Int, E::Int, type::Type = Float64
) where V <: AbstractArray{<:Number, 1}

  W         = Matrix{type}(undef, Q, E)
  P         = Matrix{Tensor{2, 3, type, 9}}(undef, Q, E)
  A         = Matrix{Tensor{4, 3, type, 81}}(undef, Q, E)
  state_old = Matrix{typeof(material_state)}(undef, Q, E)
  state_new = Matrix{typeof(material_state)}(undef, Q, E)
  # TODO move to barrier function
  for e in axes(state_old, 2)
    for q in axes(state_old, 1)
      state_old[q, e] = material_state
      state_new[q, e] = material_state
    end
  end
  S = StructArray{SectionQuadraturePointValues{length(material_state), type}}(
    (W=W, P=P, A=A, state_old=state_old, state_new=state_new)
  )
  return PreAllocatedTotalLagrangeSectionValues{length(material_state), type, typeof(S)}(S)
end

###############################################################################

# Sections
abstract type AbstractSection{FSpace, SecForm, SecVals, Mat} end
connectivity(sec::Sec) where Sec <: AbstractSection = connectivity(sec.fspace)
connectivity(sec::Sec, e::Int) where Sec <: AbstractSection = connectivity(sec.fspace, e)
dof_connectivity(sec::Sec) where Sec <: AbstractSection = dof_connectivity(sec.fspace)
dof_connectivity(sec::Sec, e::Int) where Sec <: AbstractSection= dof_connectivity(sec.fspace, e)
# shape_function_values(sec::Sec) where Sec <: AbstractSection = shape_function_values(sec.fspace)
shape_function_gradients(sec::Sec) where Sec <: AbstractSection = shape_function_gradients(sec.fspace)

struct TotalLagrangeSection{
  FSpace <: AbstractFunctionSpace, 
  SecForm <: AbstractSectionFormulation,
  SecVals <: AbstractSectionValues,
  Mat <: ConstitutiveModels.MechanicalModel,
  V <: AbstractArray{<:Number, 1}
} <: AbstractSection{FSpace, SecForm, SecVals, Mat}

  fspace::FSpace
  section_formulation::SecForm
  section_values::SecVals
  material_model::Mat
  material_props::V
end

function TotalLagrangeSection(
  coords::M,
  block::B,
  formulation::Form,
  material_model::Mat,
  material_properties::V1,
  material_state::V2
) where {M <: AbstractMatrix, B <: Block, 
         Form <: AbstractSectionFormulation, 
         Mat <: ConstitutiveModels.MechanicalModel,
         V1 <: AbstractArray,
         V2 <: AbstractArray}

  # TODO make an optional type later on
  @warn "Defaulting to second order quadrature"
  fspace = PreAllocatedFunctionSpace(coords, block, 2, 2)
  # TODO make an optional type later one
  values = PreAllocatedTotalLagrangeSectionValues(
    material_state,
    number_of_q_points(fspace), number_of_elements(fspace)
  )
  return TotalLagrangeSection(fspace, formulation, values, material_model, material_properties)
end


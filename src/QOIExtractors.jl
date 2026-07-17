abstract type AbstractQOIExtractor end

# different types to dispatch on
const H1F_QOI = 1 # H1Field
const SCM_QOI = 2 # scalar L2Field

# for component extraction
const COMPONENT_TO_INDEX = Dict{Symbol, Int}(
    :x => 1,
    :y => 2,
    :z => 3,
    # TODO figure out the tensor vs. symm tensor stuff
    # :xx => 1,
    # :yy => 2,
    # :zz => 3,
    # #
)

struct ScalarQOIExtractor{
    T,
    F <: Function,
    R <: Function
} <: AbstractQOIExtractor
    component::Union{Nothing, Int}
    func::F
    reduction::R

    function ScalarQOIExtractor(type, c, f, r)
        new{type, typeof(f), typeof(r)}(c, f, r)
    end
end

function ScalarQOIExtractor(
    func, reduction, type::Symbol;
    component::Union{Nothing, Symbol} = nothing
)
    if type == :h1_field
        # storage = create_field(objective)
        type = H1F_QOI
    elseif type == :scalar_material
        @assert component === nothing
        # storage = deepcopy(objective.assembler.scalar_quadrature_storage)
        # fill!(storage, zero(eltype(storage)))
        type = SCM_QOI
    else
        @assert false "Supported types are :h1_field and :scalar_material"
    end

    if component !== nothing
        component = COMPONENT_TO_INDEX[component]
    end

    return ScalarQOIExtractor(type, component, func, reduction)
end

function create_scratch(::ScalarQOIExtractor{H1F_QOI}, objective)
    return create_field(objective)
end

function create_scratch(::ScalarQOIExtractor{SCM_QOI}, objective)
    return deepcopy(objective.assembler.scalar_quadrature_storage)
end

function value!(val, scratch, qoi::ScalarQOIExtractor{H1F_QOI}, o, u, p)
    asm = o.assembler
    assemble_vector!(scratch, asm.vector_pattern, asm.dof, qoi.func, u, p)
    if qoi.component !== nothing
        @assert false
    else
        val[1] = reduce(qoi.reduction, scratch)
    end
    return nothing
end

function value!(val, scratch, qoi::ScalarQOIExtractor{SCM_QOI}, o, u, p)
    asm = o.assembler
    FEC.assemble_scalar_enzyme_safe!(scratch, nothing, asm.dof, qoi.func, u, p)
    val[1] = reduce(qoi.reduction, scratch)
    return nothing
end

function value(qoi::ScalarQOIExtractor, o, u, p)
    val = zeros(1)
    scratch = create_scratch(qoi, o)
    value!(val, scratch, qoi, o, u, p)
    return val[1]
end

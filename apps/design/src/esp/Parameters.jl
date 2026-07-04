struct ConstantParameters{
    T <: Number,
    V <: AbstractVector{T}
}
    names::Vector{String}
    values::V
end
Base.getindex(params::ConstantParameters, i::Int) = params.values[i]
Base.length(params::ConstantParameters) = length(params.values)

struct DesignParameters{
    T <: Number,
    V <: AbstractVector{T}
}
    names::Vector{String}
    initials::V
    lbs::V
    ubs::V
    values::V
end
Base.getindex(params::DesignParameters, i::Int) = params.values[i]
Base.length(params::DesignParameters) = length(params.values)

struct Parameters{
    T <: Number,
    V <: AbstractVector{T}
}
    constant_params::ConstantParameters{T, V}
    design_params::DesignParameters{T, V}
end

function Parameters(csm_file::String, T::Type{<:Number} = Float64)
    const_names = String[]
    const_values = T[]

    design_names = String[]
    design_values = T[]
    design_lbs = T[]
    design_ubs = T[]
    open(csm_file, "r") do f
        for line in eachline(f)
            if occursin("conpmtr", line)
                parts = split(line, " ")
                push!(const_names, parts[2])
                push!(const_values, parse(T, parts[3]))
            end

            if occursin("despmtr", line)
                parts = split(line, " ")
                push!(design_names, parts[2])
                push!(design_values, parse(T, parts[3]))
                push!(design_lbs, parse(T, parts[5]))
                push!(design_ubs, parse(T, parts[7]))
            end
        end
    end
    const_params = ConstantParameters(const_names, const_values)
    design_params = DesignParameters(design_names, design_values, design_lbs, design_ubs, design_values)

    return Parameters(const_params, design_params)
end

function Base.show(io::IO, params::Parameters; pad::String = "")
    println(io, pad, "ESP Parameters:")
    println(io, pad, "  Constant parameters:")
    for n in 1:length(params.constant_params)
        println(io, pad, "    $(params.constant_params.names[n]) = $(params.constant_params.values[n])")
    end
    println(io, pad, "  Design parameters:")
    for n in 1:length(params.design_params)
        println(io, pad, "    $(params.design_params.names[n]) = $(params.design_params.values[n])")
    end
end

function update!(params::Parameters{T, V}, design_param_vals::AbstractVector{T}) where {T, V}
    @assert length(params.design_params) == length(design_param_vals)
    for n in axes(design_param_vals, 1)
        @assert design_param_vals[n] >= params.design_params.lbs[n]
        @assert design_param_vals[n] <= params.design_params.ubs[n]
        params.design_params.values[n] = design_param_vals[n]
    end
    return nothing
end

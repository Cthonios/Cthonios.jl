struct StructuralOptimization{
    N,
    T <: Number,
    V <: AbstractVector{T},
    M <: AbstractMatrix{T},
    D <: AbstractVector{V}
} <: AbstractDesignOptimization
    f::Base.RefValue{T}
    df::V
    dfdX::H1Field{T, V, N}
    dXdp::M
    geometry::ESP.Geometry{T, V, D}
end

function StructuralOptimization(csm_file::String; esp_mode::String = "python")
    ndim = 2
    geometry = ESP.Geometry(csm_file; mode = esp_mode)
    p = copy(ESP.get_design_parameters(geometry))
    ESP.update!(geometry, p)
    p1 = reshape(geometry.dxdps[1], 3, length(geometry.dxdps[1]) ÷ 3)
    f = Base.RefValue(zero(Float64))
    df = zeros(Float64, length(p))
    dfdX = H1Field(zeros(Float64, ndim, size(p1, 2)))
    dXdp = zeros(Float64, length(p), ndim * size(p1, 2))
    return StructuralOptimization(f, df, dfdX, dXdp, geometry)
end

function get_design_parameters(opt::StructuralOptimization)
    return ESP.get_design_parameters(opt.geometry)
end

function get_design_sensitivities(opt::StructuralOptimization)

end

struct ScalarQOIExtractor{
    Assembler, 
    Func, 
    Reduction, 
    Storage
} <: AbstractQOIExtractor
    assembler::Assembler
    func::Func
    reduction::Reduction
    storage::Storage
end

function qoi_value!(f, qoi, Uu, p)
    assemble_scalar_for_ad!(qoi.storage, qoi.assembler, Uu, p, qoi.func)
    fill!(f, qoi.reduction(map(qoi.reduction, qoi.storage)))
    return nothing
end

function qoi_gradient_and_value!(f, dfdX, qoi::ScalarQOIExtractor, Uu, p)
    df = make_zero(f)
    fill!(df, one(eltype(df)))

    dqoi = make_zero(qoi)
    dUu = make_zero(Uu)
    dp = make_zero(p)

    _qoi_gradient_and_value!(f, df, qoi, dqoi, Uu, dUu, p, dp)
    copyto!(dfdX, dp.h1_coords)
    return nothing
end

function _qoi_gradient_and_value!(f, df, qoi, dqoi, Uu, dUu, p, dp)
    autodiff(
        Reverse, 
        qoi_value!,
        Duplicated(f, df),
        Duplicated(qoi, dqoi),
        Duplicated(Uu, dUu),
        Duplicated(p, dp)
    )
    return nothing
end

module CthoniosAMDGPUExt

import AMDGPU.rocSPARSE: ic0, ROCSparseMatrixCSC, ROCSparseMatrixCSR, sv2
using Adapt
using AMDGPU
using Cthonios
using KernelAbstractions
using LinearAlgebra

function Cthonios._cholesky(A::ROCSparseMatrixCSC; shift=0.0)
    A_csr = ROCSparseMatrixCSR(A)
    # TODO implement shift
    return ic0(A_csr, 'O')
end

function Cthonios._ldiv!(y, A, v, ::ROCBackend)
    temp = sv2('N', 'L', 'N', A.preconditioner, v, 'O')
    copyto!(y, temp)
    return nothing
end

end # module
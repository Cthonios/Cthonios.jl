module CthoniosAMDGPUExt

# import AMDGPU.rocBLAS: trsv!
import AMDGPU.rocSPARSE: ic0, ROCSparseMatrixCSC, ROCSparseMatrixCSR, sv2
import AMDGPU.rocSPARSE: rocsparse_create_handle
using Adapt
using AMDGPU
using Cthonios
using KernelAbstractions
using LinearAlgebra

# const handle = AMDGPU.rocSPARSE.create_handle()

function Cthonios._cholesky(A::ROCSparseMatrixCSC; shift=0.0)
    A_csr = ROCSparseMatrixCSR(A)
    # TODO implement shift
    # return ic0(A_csr, 'O')
    return ic0(A_csr, 'O')
end

function Cthonios._ldiv!(y, A, v, ::ROCBackend)
    copyto!(y, v)
    AMDGPU.rocSPARSE.sv2!('N', 'L', 'N', 1.0, A.preconditioner, y, 'O')
    AMDGPU.rocSPARSE.sv2!('T', 'L', 'N', 1.0, A.preconditioner, y, 'O')
    return nothing
end

end # module
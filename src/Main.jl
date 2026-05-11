const N = 3

function _cthonios_main(args::Vector{String})
    #####################################
    # need to define some types
    #####################################
    # N = 1 # number of fields to solve for in app
    ET = ExodusDatabase{Int32, Int32, Int32, Float64}
    SFT = AT.ScalarExpressionFunction{Float64}
    VFT = AT.VectorExpressionFunction{N, Float64}
    SPT = FEC.CSCMatrix()

    ##################################################
    # Setup app
    ##################################################
    app = AT.App{N}("Cthonios")
    sim = AT.setup(app, ARGS)

    #####################################
    # setup function space
    #####################################
    V = FunctionSpace{true}(sim.mesh, H1Field, Lagrange)
    u = ScalarFunction(V, "u")
    dof = DofManager{false}(u)
    asm = SparseMatrixAssembler{SPT, false, false}(dof)

end

function cthonios_main()::Cint
    _cthonios_main(ARGS)
    return 0
end

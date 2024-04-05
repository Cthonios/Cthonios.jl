cli_options = ArgParseSettings()
@add_arg_table! cli_options begin
  "--input-file", "-i"
    help = "an absolute path to an input file <input-file>"
    arg_type = String
    required = true
  "--ka-backend"
    help = "a backend to use with KernelAbstractions.jl such as CPU, CUDA, etc."
    arg_type = String
    default = "NoKABackend"
  "--verbose"
    action = :store_true
    help = "a flag to print to console rather than a log file"
end

# """
# $(TYPEDSIGNATURES)
# """
# function build_executable(; build_dir::String = dirname(dirname(pathof(@__MODULE__))), force::Bool = false)
#   create_app(
#     build_dir, "cthonios";
#     executables=[
#       "cthonios" => "julia_main"
#     ],
#     force=force
#     # precompile_execution_file="precompile/precompile.jl"
#   )
# end

# function dump_dependencies_state()
#   deps = Pkg.dependencies()

#   @info "Artifacts:"
#   for (uuid, dep) in deps
#     if occursin("_jll", dep.name)
#       @info "$uuid $(rpad(dep.name, 32, ' ')) $(dep.version)"
#     end
#   end
#   @info "\n"
#   @info "Dependencies:"
#   for (uuid, dep) in deps
#     if !occursin("_jll", dep.name)
#       @info "$uuid $(rpad(dep.name, 32, ' ')) $(dep.version)"
#     end
#   end
# end

function problems_main(input_file, common)
  input_settings = parse_input_file(input_file)
  prob = nothing
  for prob_settings in input_settings[:problems]
    new_section("Problem $(prob_settings[:type])")
    @timeit timer(common) "Problem $(prob_settings[:type])" begin
      type = eval(Meta.parse(prob_settings[:type]))
      prob = type(prob, prob_settings, common)
      solve!(prob, common)
    end
  end
end

"""
$(TYPEDSIGNATURES)
"""
function cthonios_main(input_file::String, verbose::Bool, ka_backend_str::String)
  log_file_name = splitext(input_file)[1] * ".log"

  # backend setup
  ka_backend = CthoniosBackend(eval(Meta.parse(ka_backend_str))())

  common = CthoniosCommon(log_file_name, ka_backend)
  cthonios_header(common)
  dump_input_file(common, input_file)

  if verbose
    # maybe add a debug flag or something like that?
    # dump_dependencies_state()
    problems_main(input_file, common)
    @info timer(common)
  else
    with_logger(common) do
      # dump_dependencies_state()
      problems_main(input_file, common)
      new_section("Timings")
      @info timer(common)
    end
  end
end

# main function, eventually add a CLI wrapper
function julia_main()::Cint
  parsed_args = parse_args(ARGS, cli_options)

  # unpack args
  input_file = parsed_args["input-file"]
  verbose = parsed_args["verbose"]
  ka_backend = parsed_args["ka-backend"]
  
  # run Cthonios
  cthonios_main(input_file, verbose, ka_backend)
  
  return 0
end

function cthonios_main_mpi end
function julia_main_mpi end

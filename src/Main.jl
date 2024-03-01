cli_options = ArgParseSettings()
@add_arg_table! cli_options begin
  "--input-file", "-i"
    help = "an absolute path to an input file <input-file>"
    arg_type = String
    required = true
  "--verbose"
    action = :store_true
    help = "a flag to print to console rather than a log file"
end

function problems_main(input_file, common)
  input_settings = parse_input_file(input_file)
  for (key, prob_settings) in input_settings[:problems]
    new_section("Problem $key")
    @timeit timer(common) "Problem $key" begin
      type = eval(Meta.parse(prob_settings[:type]))
      prob = type(prob_settings, common)
      solve!(prob, common)
    end
  end
end

function cthonios_main(input_file::String, verbose::Bool)
  log_file_name = splitext(input_file)[1] * ".log"

  common = CthoniosCommon(log_file_name)

  cthonios_header(common)
  dump_input_file(common, input_file)

  if verbose
    problems_main(input_file, common)
    @info timer(common)
  else
    with_logger(common) do
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

  # run Cthonios
  cthonios_main(input_file, verbose)
  
  return 0
end

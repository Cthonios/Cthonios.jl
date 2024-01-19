cli_options = ArgParseSettings()
@add_arg_table! cli_options begin
  "--input-file"
    help = "an absolute path to an input file <input-file>"
    arg_type = String
    required = true
end

# main function, eventually add a CLI wrapper
function julia_main()::Cint
  parsed_args = parse_args(ARGS, cli_options)
  input_file = parsed_args["input-file"]
  log_file_name = splitext(input_file)[1] * ".log"

  common = CthoniosCommon(log_file_name)

  cthonios_header(common)
  dump_input_file(common, input_file)

  with_logger(common) do
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

  with_logger(common) do
    new_section("Timings")
    @info timer(common)
  end

  return 0
end

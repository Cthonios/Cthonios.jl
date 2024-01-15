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

  # input_settings = YAML.load_file(input_file)
  # # domain = Domain(common, input_file)

  # # domain = setup(QuasiStaticDomain, common, input_settings)
  # @assert "problem" in keys(input_settings)
  # @assert "type" in keys(input_settings["problem"])
  # problem_type = Meta.parse(input_settings["problem"]["type"])
  # with_logger(common) do
  #   new_section("Problem setup")
  #   problem = eval(problem_type)(input_settings)
  #   new_section("Problem solution")
  #   solve!(problem)
  # end

  with_logger(common) do
    input_settings = parse_input_file(input_file)
    for prob_settings in values(input_settings[:problems])
      type = eval(Meta.parse(prob_settings[:type]))
      prob = type(prob_settings)
      solve!(prob)
    end
  end

  # new_section(common, "Timings")
  with_logger(common) do
    new_section("Timings")
    @info timer(common)
  end

  return 0
end


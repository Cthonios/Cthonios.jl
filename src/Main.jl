function cthonios_main()::Cint
  settings = parse_command_line()
  input_settings = read_input_file(settings)
  run(input_settings)
  return Int32(0)
end

function parse_command_line()::Dict{String, Any}
  settings = ArgParseSettings()
  @add_arg_table! settings begin
    "--input-file", "-i"
      help = "Path to an input file"
      arg_type = String
  end
  return parse_args(settings)
end

function read_input_file(settings)
  input_file = settings["input-file"]
  YAML.load_file(input_file; dicttype=Dict{Symbol, Any})
end

function run(input_settings)
  problems = input_settings[:problems]
  for prob_inputs in problems
    type = Symbol(prob_inputs[:type])
    problem, p = @eval $type($prob_inputs)
    Uu = create_unknowns(problem.solver)
    solve!(problem, Uu, p)
  end
end

function cthonios_main()::Cint
  cli_args = parse_command_line()
  @assert cli_args["backend"] in ["cpu", "cuda", "rocm"]
  @assert cli_args["input-file"] !== nothing
  @info "cthonios_main"
  @info "Backend    = $(cli_args["backend"])"
  @info "Input file = $(cli_args["input-file"])"
  input_settings = parse_input_file(cli_args)
  @show input_settings
  return Int32(0)
end

function parse_command_line()::Dict{String, Any}
  settings = ArgParseSettings()
  @add_arg_table! settings begin
    "--backend", "-b"
      arg_type = String
      default = "cpu"
      help = "Backend to use, e.g. cpu, cuda, rocm"
    "--input-file", "-i"
      arg_type = String
      help = "Path to an input file"
      required = true
  end
  return parse_args(settings)
end

function parse_input_file(settings)
  input_file = settings["input-file"]
  @info "parsing input file = $(input_file)"
  YAML.load_file(input_file; dicttype=Dict{Symbol, Any})
end

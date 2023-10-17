# some exceptions for cli error handling
# function cthonios_cli_error(io::IO, error_message::String)
#   print(io, "\n\nCthonios error:\n\n", error_message, "\n")
# end

struct InputFileNotFoundException <: CthoniosException
  parsed_args::Dict{String, Any}
end

function Base.show(io::IO, e::InputFileNotFoundException)
  str = "Input file not found error!\n"
  str = str * "\nInput file: \"$(e.parsed_args["input-file"])\" "
  str = str * "not found in directory \"$(abspath(dirname(e.parsed_args["input-file"])))\"\n"
  # cthonios_cli_error(io, str)
  print(io, "\n\nCthonios error:\n\n", error_message, "\n")
end

function input_file_not_found(logger::SimpleLogger, parsed_args::Dict{String, Any})
  e = InputFileNotFoundException(parsed_args)
  with_logger(logger) do 
    str_e = string(e)
    @error str_e
  end
  throw(e)
end

# main function for compiled code
function internal_main(parsed_args::Dict{String, <:Any})
  
  # start logging
  log_file_name = splitext(abspath(parsed_args["input-file"]))[1] * ".log"
  logger = FormatLogger(open(log_file_name, "w")) do io, args
    println(io, args.message)
  end

  # print the header and start doing some error checking on command line arguments
  with_logger(logger) do 
    @info cthonios_header(parsed_args)

    input_file = aprepro(parsed_args["input-file"])
    if !isfile(abspath(input_file))
      input_file_not_found(logger, parsed_args)
    end

    input_file_lines = readlines(input_file)
    input_file_str = section_header("Input file")
    input_file_str = input_file_str * "Input file = $input_file\n\n" 
    for line in input_file_lines
      input_file_str = input_file_str * line * "\n"
    end
    input_file_str = input_file_str * "\n\n"
    @info input_file_str

    input_settings = input_file_parser(input_file)

    System(input_settings)

  end

  # cleanup
  # close(log_file_io)
  # return Cint(0)
end

internal_main(file_name::String) = 
internal_main(Dict("input-file" => file_name))

function julia_main()::Cint
  parsed_args = cli_parser()
  internal_main(parsed_args)
  return Cint(0)
end
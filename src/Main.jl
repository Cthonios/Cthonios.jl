# main function, eventually add a CLI wrapper
function cthonios_main(input_file::String)
  log_file_name = splitext(input_file)[1] * ".log"

  common = CthoniosCommon(log_file_name)

  cthonios_header(common)
  dump_input_file(common, input_file)

  input_settings = YAML.load_file(input_file)
  # domain = Domain(common, input_file)

  domain = setup(QuasiStaticDomain, common, input_settings)

  # new_section(common, "Timings")
  with_logger(common) do
    new_section("Timings")
    @info timer(common)
  end

  return domain
end
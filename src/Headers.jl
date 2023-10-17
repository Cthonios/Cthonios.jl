# some printing stuff
function section_header(section_header::String)
  str = repeat("=", 80) * "\n"
  str = str * "= $section_header \n"
  str = str * repeat("=", 80) * "\n\n"
  return str
end

function cthonios_header(parsed_args::Dict{String, <:Any})
  str = raw"""

   ____     __    __                                                   ___      
  /\  _`\  /\ \__/\ \                      __                      __ /\_ \     
  \ \ \/\_\\ \ ,_\ \ \___     ___     ___ /\_\    ___     ____    /\_\\//\ \    
   \ \ \/_/_\ \ \/\ \  _ `\  / __`\ /' _ `\/\ \  / __`\  /',__\   \/\ \ \ \ \   
    \ \ \L\ \\ \ \_\ \ \ \ \/\ \L\ \/\ \/\ \ \ \/\ \L\ \/\__, `\__ \ \ \ \_\ \_ 
     \ \____/ \ \__\\ \_\ \_\ \____/\ \_\ \_\ \_\ \____/\/\____/\_\_\ \ \/\____\
      \/___/   \/__/ \/_/\/_/\/___/  \/_/\/_/\/_/\/___/  \/___/\/_/\ \_\ \/____/
                                                                  \ \____/      
                                                                   \/___/       
  """
  str = str * "\n\nCthonios version: α0.1.0\n\n"
  str = str * "Continuum Thermodynamics Optimization-based\n"
  str = str * "Numerical Implementations of Solid mechanics\n\n\n"
  str = str * "Developed by Craig M. Hamel\n\n\n"
  str = str * section_header("Command line arguments")
  max_key_length = length.(keys(parsed_args)) |> maximum
  for (key, val) in parsed_args
    str = str *"  $(rpad(key, max_key_length)) => $val\n"
  end
  str = str * "\n\n"
  return str
end

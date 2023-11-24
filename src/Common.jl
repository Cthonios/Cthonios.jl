################################################### 
# Logger below, for writing to log file
################################################### 

# TODO overhaul logger to have a notion of loggin levels
# to automagically handle indentation for you
"""
"""
struct CthoniosLogger
  logger::FormatLogger
end

"""
"""
function CthoniosLogger(file_name::String)
  logger = FormatLogger(file_name) do io, args
    if args.level == Error
      println(io, "[", args.level, "] ", args.message)
    elseif args.level == Warn
      println(io, "[", args.level, "] ", args.message)
    else
      println(io, args.message)
    end
  end
  return CthoniosLogger(logger)
end

"""
"""
function with_logger(f, logger::CthoniosLogger)
  Base.with_logger(f, logger.logger)
end



################################################### 
# Timer below
################################################### 

"""
"""
struct CthoniosTimer
  timer::TimerOutput
end

"""
"""
function CthoniosTimer()
  return CthoniosTimer(TimerOutput())
end

Base.show(io::IO, timer::CthoniosTimer) = Base.show(io::IO, timer.timer)

################################################### 
# Bundled struct of common stuff
################################################### 

# TODO I'm sure more will be added, like post-processing
# or parsing, etc.

"""
"""
struct CthoniosCommon
  logger::CthoniosLogger
  timer::CthoniosTimer
end

"""
"""
function CthoniosCommon(file_name::String)
  return CthoniosCommon(
    CthoniosLogger(file_name), 
    CthoniosTimer()
  )
end

"""
"""
function with_logger(f, common::CthoniosCommon)
  with_logger(f, common.logger)
end

# function with_timer(f, common::CthoniosCommon, section::String)
#   @timeit common.timer.timer section begin
#     f
#   end
# end

timer(common::CthoniosCommon) = common.timer.timer

# new_section(common::CthoniosCommon, section_name::String) = 
# new_section(common.logger, section_name)

"""
"""
function new_section(section_name::String)
  string = "\n" * repeat('=', 80)
  string = string * "\n= $section_name\n"
  string = string * repeat('=', 80)
  # with_logger(logger) do
  #   @info string
  # end
  @info string
end

function cthonios_header(common::CthoniosCommon)
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
  
  with_logger(common) do
    @info str
  end
end

function dump_input_file(common::CthoniosCommon, file_name::String)
  lines = readlines(file_name)

  with_logger(common) do 
    @info "Reading from input file $(abspath(file_name))\n"
    new_section("Input FIle")
    for line in lines
      @info line
    end
    # @info "\n"
  end
end

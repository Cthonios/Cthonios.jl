# struct VariableTable
#   variables::
# end

struct VariableList
  variables::Vector{Symbol}
end

VariableList(variable_names::Vector{String}) = VariableList(Symbol.(variable_names))

# make a list of variables
# add variabl name to kernel input file stuff
# add coupling variable names to kernel input file stuff
# finally make a simple binary table showing couplings
# TODO
# plans for QOI extractors
# 1. have a base cache type to handle whether
# or not we need AD which is the main
# source of allocations for these structs
# 2. have input helper methods for what type of 
# calc is needed, how to reduce, and what types
# to pre-allocate for storage

abstract type AbstractQOIExtractor end

function value end
function value! end

include("ScalarQOIExtractor.jl")

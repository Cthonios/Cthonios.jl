# struct IteratorType{T}
#   type::T
# end
abstract type AbstractIteratorType end

struct AKIteator
end

struct VanillaIterator
end

include("InPlaceIterators.jl")
include("OutOfPlaceIterators.jl")

# hopefully below will be default soon
# include("AKIterators.jl")

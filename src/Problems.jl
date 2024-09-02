abstract type AbstractProblem end

struct Problem{S}
  solver::S
end

function solve!(prob::Problem)

end

# exports
export Problem
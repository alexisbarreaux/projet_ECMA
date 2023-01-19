using JuMP
using CPLEX

include("../../utils/constants.jl")
include("../../utils/jsonUtils.jl")
include("../../utils/instancesUtils.jl")

"""
include("models/subModels/partWeightsModel.jl")
computeAndStoreBoundAllFiles("weightBounds2.json")
"""
function computeBoundAllFiles()::Dict{String, Float64}
    bounds = Dict{String, Float64}()
    for inputFile in DATA_FILES
        bounds[inputFile] = round(computeLowerBound(inputFile))
    end
    return bounds
end

function computeAndStoreBoundAllFiles(resultFile::String="test.json")::Nothing
    bounds = computeBoundAllFiles()
    filePath =RESULTS_DIR_PATH * "\\" * resultFile
    jsonDropToFile(filePath, bounds)
end

function computeLowerBound(inputFile::String)::Float64
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    gap = sum(w_v) - (K-1) * B
    return max(0, gap)
end


function partWeightsSolve(inputFile::String)::Any
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    l = computeDistances(coordinates)

    # Creating the model
    model = Model(CPLEX.Optimizer)
    set_silent(model)

    # Variables
    @variable(model, y[i in 1:n, k in 1:K], Bin)

    # Objective
    @objective(model, Min, sum(y[i,1] for i in 1:n))
    
    # Constraints
    # Weights of the parts
    @constraint(model, [k in 1:K], sum(w_v[i]*y[i,k] for i in 1:n) <= B)

    # Each node is in a part
    @constraint(model, [i in 1:n],  sum(y[i,k] for k in 1:K) == 1)

    # We can also add that, because of the base problem, each part must have at least one node
    @constraint(model, [k in 1:K],  sum(y[i,k] for i in 1:n) >= 1)

    # Solve
    optimize!(model)

    ### Display the solution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
    # Récupération des valeurs d’une variable
        result = JuMP.value.(y)
        value = JuMP.objective_value(model)
        #println("Success, nodes : " * string(JuMP.node_count(model))* ", Time : "* string(round(JuMP.solve_time(model), digits= 5)) * " Value : " * string(round(value, digits=4)))
        for k in 1:K
            partWeight = 0
            for i in 1:n
                if result[i,k] == 1
                    partWeight += w_v[i]
                end
            end
            #println("Part " * string(k) * " weight " * string(partWeight))
        end
        #println("B " * string(B))
        println(inputFile, " ", value)
        return value
    else
        println("Not feasible!!")
        return
    end

end

computeAndStoreBoundAllFiles()
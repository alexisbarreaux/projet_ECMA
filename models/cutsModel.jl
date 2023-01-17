using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/instancesUtils.jl")
include("../utils/jsonUtils.jl")

# include("models/cutsModel.jl")
# cutsSolve("10_ulysses_3.tsp")
function solveAndReturnAllInstancesCuts()::Dict{String, Float64}
    cutsValues = Dict{String, Float64}()
    for inputFile in DATA_FILES
        cutsValues[inputFile] = cutsSolve(inputFile)
    end
    return cutsValues
end

function solveAndStoreAllInstancesStatic(resultFile::String=DUAL_RESULTS_FILE)::Nothing
    cutsValues = solveAndReturnAllInstancesCuts()
    filePath =RESULTS_DIR_PATH * "\\" * resultFile
    jsonDropToFile(filePath, cutsValues)
end


function cutsSolve(inputFile::String, showResult::Bool= false, silent::Bool=true)::Any
    """
    The source file includes the following variables:
        - n : number of nodes,
        - L : bound on the total uncertainty on the lengths,
        - W : bound on the total uncertainty for the weights,
        - K : number of parts in which to split the graph,
        - B : capacity constraints on the parts,
        - w_v : weights of the nodes,
        - W_v : bound on the uncertainty of the weight for each node
        - lh : lengths linked to the uncertainty on lengths for each node.
        - coordinates : coordinates of our points/nodes.
    """
    println("Solving ", inputFile, " in cuts mode.")
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    l = computeDistances(coordinates)

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    # Variables
    @variable(model, x[i in 1:n, j in 1:n], Bin)
    @variable(model, y[i in 1:n, k in 1:K], Bin)
    # New variables for cuts
    @variable(model, z >= 0)
    # Also "variables" but not handled by the model
    U_1 = [l]
    U_2 = [w_v]

    # Objective
    @objective(model, Min, z)
    
    # Triangular inequalities between x and y
    @constraint(model, [k in 1:K, i in 1:n, j in i:n], y[i,k] + y[j,k] <= x[i,j] + 1)
    # Each node is in a part
    @constraint(model, [i in 1:n],  sum(y[i,k] for k in 1:K) == 1)
    # New constraints for cut
    # Rewriting of the objective
    @constraint(model, [l_1 in U_1],z >=  sum(x[i,j] * l_1[i,j]))
    # Weights of the parts
    @constraint(model, [w_2 in U_2], sum(w_2[i]*y[i,k] for i in 1:n) <= B)

    # Solve
    start = time()
    optimize!(model)
    optimize_time = time() - start

    ### Display the solution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
    # Récupération des valeurs d’une variable
        result = JuMP.value.(y)
        value = JuMP.objective_value(model)
        
        if showResult
            println("Success, nodes : " * string(JuMP.node_count(model))* ", Time : "* string(round(optimize_time, digits= 5)) * " Value : " * string(round(value, digits=4)))
            createdParts = Dict{Int, Array}(k => [] for k in 1:K)
            for i in 1:n
                for k in 1:K
                    if result[i,k] == 1
                        createdParts[k] = vcat(createdParts[k],[i])
                        break
                    end
                end
            end
            println("Found parts are : ", createdParts)
        end
        return value
    else
        println("Not feasible!!")
        return
    end

end
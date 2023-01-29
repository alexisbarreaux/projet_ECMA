using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/instancesUtils.jl")

"""
include("models/staticModel.jl")
staticSolve("10_ulysses_3.tsp")
"""
# Expected result :  {1, 5, 8} {2, 3, 4}, {6, 7, 9, 10} 
# Expected value : 54.354823588

function staticSolve(inputFile::String, showResult::Bool= false, silent::Bool=true, timeLimit::Float64=60.0)::Any
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
    println("Solving ", inputFile, " in static mode.")
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    l = computeDistances(coordinates)

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if timeLimit >= 0
        set_time_limit_sec(model, timeLimit)
    end

    if silent
        set_silent(model)
    end

    # Variables
    @variable(model, x[i in 1:n, j in i+1:n], Bin)
    @variable(model, y[i in 1:n, k in 1:K], Bin)

    # Objective
    @objective(model, Min, sum(x[i,j] * l[i,j] for i in 1:n for j in i+1:n))
    
    # Constraints
    # Weights of the parts
    @constraint(model, [k in 1:K], sum(w_v[i]*y[i,k] for i in 1:n) <= B)
    # Triangular inequalities between x and y
    @constraint(model, [k in 1:K, i in 1:n, j in i+1:n], - x[i,j] + y[i,k] + y[j,k] <=1)
    # Note : adding the others constraints below doesn't seem to yield a better results in times or nodes.
    #@constraint(model, [k in 1:K, i in 1:n, j in i+1:n], x[i,j] + y[i,k] - y[j,k] <=1)
    #@constraint(model, [k in 1:K, i in 1:n, j in i+1:n], x[i,j] - y[i,k] - y[j,k] <=1)

    # Each node is in a part
    @constraint(model, [i in 1:n],  sum(y[i,k] for k in 1:K) == 1)

    # First node can be put anywhere
    # Note : this seem to actually slow down the process
    #@constraint(model, y[1,1] ==1)

    # Each part can't contain more than maxPartSize elements
    # Note : this seem to add a slight overhead and no real gain.
    #maxPartSize = boundOnPartsSize(B, w_v)
    #@constraint(model, [k in 1:K],  sum(y[i,k] for i in 1:n) <= maxPartSize)

    # Solve
    optimize!(model)

    ### Display the solution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
    # Récupération des valeurs d’une variable
        result = JuMP.value.(y)
        value = JuMP.objective_value(model)
        solveTime = round(JuMP.solve_time(model), digits= 5)
        gap = JuMP.relative_gap(model)
        if showResult
            println("Success, nodes : " * string(JuMP.node_count(model))* ", Time : "* string(solveTime) * " Value : " * string(round(value, digits=4)))
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
        return isOptimal, solveTime, value, gap
    else
        println("Not feasible!!")
        return
    end

end

function staticSolveUnoptimizedInstance(timeLimit::Float64=60.0, resultFile::String=STATIC_RESULTS_FILE)::Nothing
    # Loading
    filePath =RESULTS_DIR_PATH * "\\" * resultFile * ".csv"
    # Get unoptimal instance
    if !isfile(filePath)
        println("No such file to load")
        return
    else
        currentResults = DataFrame(CSV.File(filePath))
        fileToRun, rowToReplace = findUnoptimzedInstance(currentResults)
    end

    if fileToRun == ""
        println("No non optimized instance is left.")
        return
    else
        # Run
        println("Found unoptimized instance " * fileToRun)
        if runInstanceAndUpdateDataframe(currentResults, fileToRun, timeLimit, rowToReplace)
            CSV.write(filePath, currentResults, delim=";")
        end
        return 
    end
end

function staticSolveChosenInstance(fileToRun::String, timeLimit::Float64=60.0, resultFile::String=STATIC_RESULTS_FILE)::Nothing
    # Loading
    filePath =RESULTS_DIR_PATH * "\\" * resultFile * ".csv"
    # Get unoptimal instance
    if !isfile(filePath)
        println("No such file to load")
        return
    else
        currentResults = DataFrame(CSV.File(filePath))
    end

    # Run
    if runInstanceAndUpdateDataframe(currentResults, fileToRun, timeLimit)
        CSV.write(filePath, currentResults)
    end

    return 
end

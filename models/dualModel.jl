using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/instancesUtils.jl")
include("../utils/jsonUtils.jl")

include("../heuristic/heuristic.jl")

"""
include("models/dualModel.jl")
dualSolve("10_ulysses_3.tsp")

Solution robuste : {1, 2, 3, 10}, {4, 6, 7, 8}, {5, 9}
(objectif : 136.995276296)
"""


function dualSolve(inputFile::String, showResult::Bool= false, silent::Bool=true, timeLimit::Float64=60.0, warm::Bool=false)::Any
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
    println("Solving ", inputFile, " in dual mode.")
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    l = computeDistances(coordinates)

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end
    if timeLimit >= 0
        set_time_limit_sec(model, timeLimit)
    end
    # CPLEX params
    """
    set_optimizer_attribute(model, "CPX_PARAM_CLIQUES", 3)
    set_optimizer_attribute(model, "CPX_PARAM_COVERS", 3)
    set_optimizer_attribute(model, "CPX_PARAM_ZEROHALFCUTS", 1)
    set_optimizer_attribute(model, "CPX_PARAM_MIRCUTS", 2)
    """

    # Variables
    if warm
        solution = run_first_heuristic(inputFile, true)
        if solution != nothing
            # Using heuristic
            x_heur, y_heur = translate_output(inputFile, solution)
            #_, value = heuristicResult
            #set_optimizer_attribute(model, "CPX_PARAM_CUTUP", value)
            @variable(model, x[i=1:n, j=i+1:n], Bin, start=x_heur[i,j])
            @variable(model, y[i=1:n, k=1:K], Bin, start = y_heur[i,k])
        else
            @variable(model, x[i in 1:n, j in i+1:n], Bin)
            @variable(model, y[i in 1:n, k in 1:K], Bin)
        end
    else
        @variable(model, x[i in 1:n, j in i+1:n], Bin)
        @variable(model, y[i in 1:n, k in 1:K], Bin)
    end
    
    
    # New constraints added for the dual
    @variable(model, beta[i in 1:n, j in i+1:n] >= 0.)
    @variable(model, alpha >= 0.)
    @variable(model, dzeta[i in 1:n, k in 1:K] >= 0.)
    @variable(model, gamma[k in 1:K] >= 0.)

    # Objective
    @objective(model, Min, sum((x[i,j] * l[i,j]) + 3*beta[i,j] for i in 1:n for j in i+1:n) + L*alpha)
    
    # Constraints
    # Triangular inequalities between x and y
    @constraint(model, [k in 1:K, i in 1:n, j in i+1:n], y[i,k] + y[j,k] <= x[i,j] + 1)
    # Each node is in a part
    @constraint(model, [i in 1:n],  sum(y[i,k] for k in 1:K) == 1)
    # Link between alpha, beta and x (comes from first dual on objective)
    @constraint(model, [i in 1:n, j in i+1:n], alpha + beta[i,j] >= (lh[i] + lh[j])*x[i,j])
    # Constraint on the parts from the robust one dualized
    @constraint(model, [k in 1:K], (sum(w_v[i]*y[i,k] + W_v[i] * dzeta[i,k] for i in 1:n) + W*gamma[k]) <= B)
    # Constraint from the second dual
    @constraint(model, [i in 1:n, k in 1:K], (gamma[k] + dzeta[i,k]) >= w_v[i]*y[i,k])    

    # First node can be put anywhere
    @constraint(model, y[1,1] ==1)

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
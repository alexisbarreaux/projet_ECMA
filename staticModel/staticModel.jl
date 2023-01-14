using JuMP
using CPLEX
include("../utils/constants.jl")
include("../utils/utils.jl")

# staticSolve("10_ulysses_3.tsp")

function staticSolve(inputFile::String)::String
    """
    The source file includes the following variables:
        - n : number of nodes,
        - L :
        - W :
        - K : number of parts in which to split the graph,
        - B : capacity constraints on the parts,
        - w_v : weights of the nodes,
        - W_v : 
        - lh :
        - coordinates : coordinates of our points/nodes.
    """
    println("Used input file is : ", inputFile)
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    l = computeDistances(coordinates)

    # Creating the model
    model = Model(CPLEX.Optimizer)

    # Variables
    @variable(model, x[i in 1:n, j in 1:n], Bin)
    @variable(model, y[i in 1:n, k in 1:K], Bin)

    # Objective
    # We only sum on i < j as per how the subject specified it on an example
    @objective(model, Max, sum(x[i,j] * l[i,j] for i in 1:n for j in i+1:n))
    
    # Constraints
    # Weights of the parts
    @constraint(model, [k in 1:K], sum(w_v[i]*y[i,k] for i in 1:n) <= B)
    # Triangular inequalities between x and y
    @constraint(model, [k in 1:K, i in 1:n, j in i+1:n], y[i,k] + y[j,k] <= x[i,j] + 1)
    # Each node is in a part
    @constraint(model, [i in 1:n],  sum(y[i,k] for k in 1:K) == 1)

    # Solve
    start = time()
    optimize!(model)
    optimize_time = time() - start

    ### Display the solution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
    # Récupération des valeurs d’une variable
        return "Success, nodes : " * string(JuMP.node_count(model))* ", Time : "* string(optimize_time) * " value " * string(round(Int, JuMP.objective_value(model)))
    else
        return "Not feasible"
    end
end
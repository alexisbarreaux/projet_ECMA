using JuMP
using CPLEX
include("../utils/constants.jl")

# staticSolve("10_ulysses_3.tsp")

function staticSolve(inputFile::String)::Any
    """
    The source file includes the following variables:
        - n : number of nodes,
        - L :
        - W :
        - K : number of parts in which to split the graph.
        - B :
        - w_v :
        - W_v : 
        - lh :
        - coordinates : .
    """
    println("Used input file is : ", inputFile)
    # Directly load data file
    include(DATA_DIR_PATH * "\\" * inputFile)

    # Creating the model
    model = Model(CPLEX.Optimizer)

    # Variables
    @variable(model, x[i in 1:n, j in 1:n], Bin)
    @variable(model, y[i in 1:n, k in 1:K], Bin)

    # Objective
    @objective(model, Max, sum(x[i,j] * similarities[i,j] for i in 1:n for j in i:n))
    
    
    # Contraintes
    
    return
end
using JuMP
using CPLEX

include("../utils/constants.jl")
include("../utils/instancesUtils.jl")
include("../utils/jsonUtils.jl")

"""
include("models/cutsModel.jl")
cutSolve("10_ulysses_3.tsp")

Solution robuste : {1, 2, 3, 10}, {4, 6, 7, 8}, {5, 9}
(objectif : 136.995276296)
"""

function firstSubProblem(x_val::Matrix{Float64},n::Int64, l::Matrix{Float64}, lh::Vector{Int64}, L::Int64)::Tuple{Float64, Matrix{Float64}, Float64}
    """
    Function to solve the first sub problem in the current state.
    """
    # Creating the model
    sub_model = Model(CPLEX.Optimizer)
    set_silent(sub_model)

    ##### Variables #####
    @variable(sub_model, 0. <= delta_1[i in 1:n, j in 1:n] <= 3.)

    ##### Objective #####
    @objective(sub_model, Max, sum( (l[i,j] + delta_1[i,j]* (lh[i] + lh[j]) ) * x_val[i,j] for i in 1:n for j in i+1:n))
    
    ##### Constraints #####
    @constraint(sub_model, sum( delta_1[i,j] for i in 1:n for j in i+1:n) <= L)
    
    optimize!(sub_model)

    return JuMP.objective_value(sub_model), JuMP.value.(delta_1), JuMP.solve_time(sub_model)
end

function secondKSubProblem(y_val::Matrix{Float64}, k::Int64, n::Int64, w_v::Vector{Int64}, W_v::Vector{Float64}, W::Int64)::Tuple{Float64, Vector{Float64}, Float64}
    """
    Function to solve the k-th second sub problem in the current state.
    """
    # Creating the model
    sub_model = Model(CPLEX.Optimizer)
    set_silent(sub_model)

    ##### Variables #####
    @variable(sub_model, 0. <= delta_2[i in 1:n] <= W_v[i])

    ##### Objective #####
    @objective(sub_model, Max, sum(w_v[i]* (1 + delta_2[i]) * y_val[i,k] for i in 1:n))
    
    ##### Constraints #####
    @constraint(sub_model, sum( delta_2[i] for i in 1:n) <= W)
    
    optimize!(sub_model)

    return JuMP.objective_value(sub_model), JuMP.value.(delta_2), JuMP.solve_time(sub_model)
end

function cutSolve(inputFile::String, showResult::Bool= false, silent::Bool=true, timeLimit::Float64=60.0)::Any
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

    ##### Variables #####
    @variable(model, x[i in 1:n, j in 1:n], Bin)
    @variable(model, y[i in 1:n, k in 1:K], Bin)
    # New variables for cuts
    @variable(model, z >= 0)

    ##### Objective #####
    @objective(model, Min, z)
    
    ##### Constraints #####
    # Triangular inequalities between x and y
    @constraint(model, [k in 1:K, i in 1:n, j in i+1:n], y[i,k] + y[j,k] <= x[i,j] + 1)
    # Each node is in a part
    @constraint(model, [i in 1:n], sum(y[i,k] for k in 1:K) == 1)
    # New constraints for cut
    # Rewriting of the objective
    @constraint(model,z >=  sum(x[i,j] * l[i,j] for i in 1:n for j in i+1:n))
    # Weights of the parts
    @constraint(model, [k in 1:K], sum(w_v[i]*y[i,k] for i in 1:n) <= B)

    # First node can be put anywhere
    @constraint(model, y[1,1] ==1)

    hasAddedConstraint = true
    runTime = 0
    while hasAddedConstraint && (timeLimit < 0 || runTime < timeLimit)
        hasAddedConstraint = false
        # Solve current state
        if timeLimit >= 0
            set_time_limit_sec(model, timeLimit - runTime)
        end

        optimize!(model)
        feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
        isOptimal = termination_status(model) == MOI.OPTIMAL
        runTime += JuMP.solve_time(model)

        if timeLimit >= 0 && runTime > timeLimit
            break
        end

        if feasibleSolutionFound
            value = JuMP.objective_value(model)
            # Solve sub problems with current optimum
            z_val = JuMP.value(z)
            x_val = JuMP.value.(x)
            y_val = JuMP.value.(y)

            z_1, delta_1, time_1= firstSubProblem(x_val, n, l, lh, L)
            runTime += time_1
            if z_1 > (z_val + 1e-5)
                if timeLimit >= 0 && runTime > timeLimit
                    println("Not feasible!!")
                    return
                else
                    cstr = @constraint(model, z >= sum(x[i,j] * (l[i,j] + delta_1[i,j]* (lh[i] + lh[j])) for i in 1:n for j in i+1:n))
                    hasAddedConstraint = true
                end
            end
    
            for k in 1:K 
                z_2_k, delta_2_k, time_2_k = secondKSubProblem(y_val, k, n, w_v, W_v, W)
                runTime += time_2_k
                if z_2_k > (B + 1e-5)
                    if timeLimit >= 0 && runTime > timeLimit
                        println("Not feasible!!")
                        return
                    else
                        # Adding constraints for all k because it gains time
                        cstr = @constraint(model, [k_2 in 1:K], sum( w_v[i]* (1 + delta_2_k[i]) * y[i,k_2] for i in 1:n) <= B)
                        hasAddedConstraint = true
                    end
                end
            end
        else
            break
        end

    end

    ### Display the solution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    gap = JuMP.relative_gap(model)
    if feasibleSolutionFound
    # Récupération des valeurs d’une variable
        # Il faut checker si on a encore des plans non valides
        value = JuMP.objective_value(model)
        # Solve sub problems with current optimum
        z_val = JuMP.value(z)
        x_val = JuMP.value.(x)
        y_val = JuMP.value.(y)

        z_1, delta_1= firstSubProblem(x_val, n, l, lh, L)
        if z_1 > (z_val + 1e-5)
            println("Not feasible!!")
            return
        end

        for k in 1:K 
            z_2_k, delta_2_k = secondKSubProblem(y_val, k, n, w_v, W_v, W)
            if z_2_k > (B + 1e-5)
                println("Not feasible!!")
                return
            end
        end


        result = JuMP.value.(y)
        solveTime = round(runTime, digits= 5)
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
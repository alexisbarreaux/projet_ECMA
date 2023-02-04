using DataStructures
using Random

include("../utils/constants.jl")
include("../utils/instancesUtils.jl")
include("../utils/jsonUtils.jl")

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

"""
include("heuristic/heuristic.jl")
run_heuristic("52_berlin_9.tsp")
"""

function compute_worst_weight(inputFile::String, cluster::Array)::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    combined_dict = OrderedDict(zip(cluster, [-w_v[k] for k in cluster]))
    total_weight = 0
    remaining_uncertainty = W
    for (node, _) in combined_dict
        error = min(remaining_uncertainty, W_v[node])
        remaining_uncertainty -= error
        total_weight += w_v[node] * (1 + error)
    end
    return total_weight
end

function compute_worst_weight_sort(inputFile::String, cluster::Array)::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    combined_dict = OrderedDict(zip(cluster, [-w_v[k] for k in cluster]))
    sorted_dict = sort(combined_dict; byvalue=true)
    total_weight = 0
    remaining_uncertainty = W
    for (node, _) in sorted_dict
        error = min(remaining_uncertainty, W_v[node])
        remaining_uncertainty -= error
        total_weight += w_v[node] * (1 + error)
    end
    return total_weight
end

function check_feasibility(inputFile::String, solution::Dict{Int,Array})::Bool
    include(DATA_DIR_PATH * "\\" * inputFile)
    for (_, cluster) in solution
        worst_weight = compute_worst_weight(inputFile, cluster)
        if worst_weight > B
            #println("Worst weight : ", worst_weight, ", B : ", B)
            return false
        end
    end
    return true
end

function compute_worst_case(inputFile::String, solution::Dict{Int,Array}, mode::Int=0)::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    l = computeDistances(coordinates)
    total_distance = 0
    supplementary_distances = []
    for (_, cluster) in solution
        for a in cluster
            for b in cluster
                if a < b
                    append!(supplementary_distances, lh[a] + lh[b])
                    total_distance += l[a, b]
                end
            end
        end
    end
    supplementary_distances = reverse!(sort(supplementary_distances))
    remaining_uncertainty = L
    robust_cost = 0
    for d in supplementary_distances
        error = min(3, remaining_uncertainty)
        remaining_uncertainty -= error
        robust_cost += error * d
    end
    if mode > 0
        return total_distance + robust_cost, supplementary_distances, robust_cost, l
    end
    return total_distance + robust_cost
end

function compute_score_cluster(inputFile::String, cluster::Array)::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    l = computeDistances(coordinates)
    total_distance = 0
    supplementary_distances = []
    for a in cluster
        for b in cluster
            if a != b
                append!(supplementary_distances, lh[a] + lh[b])
                total_distance += l[a, b]
            end
        end
    end
    supplementary_distances = reverse!(sort(supplementary_distances))
    remaining_uncertainty = L
    for d in supplementary_distances
        error = min(3, remaining_uncertainty)
        remaining_uncertainty -= error
        total_distance += error * d
    end
    return total_distance
end


function construct_solution(inputFile::String, mode::Int=0)::Dict{Int,Array}
    include(DATA_DIR_PATH * "\\" * inputFile)
    solution = Dict{Int,Array}(k => [] for k in 1:K)
    if mode == 0
        combined_dict = OrderedDict(zip(1:n, [-w_v[k] for k in 1:n]))
    end
    if mode == 1
        combined_dict = OrderedDict(zip(1:n, [(-w_v[k], -W_v[k]) for k in 1:n]))
    end
    if mode == 2
        combined_dict = OrderedDict(zip(1:n, [(-w_v[k], W_v[k]) for k in 1:n]))
    end
    """if mode == 3
        combined_dict = OrderedDict(zip(1:n, [-w_v[n+1-k] for k in 1:n]))
    end"""
    if mode == 3
        #Random.seed!(2)
        random_values = rand(Float64, n)
        combined_dict = OrderedDict(zip(1:n, random_values))
    end
    sorted_dict = sort(combined_dict; byvalue=true)
    for (node, _) in sorted_dict
        min_weight = 2 * B
        min_index = 0
        for k in 1:K
            append!(solution[k], node)
            if mode == 3
                worst_weight = compute_worst_weight_sort(inputFile, solution[k])
            else
                worst_weight = compute_worst_weight(inputFile, solution[k])
            end
            pop!(solution[k])
            if worst_weight < min_weight
                min_weight = worst_weight
                min_index = k
            end
        end
        if min_index == 0
            println("Heuristic solution for ", inputFile, " is not feasible in the worst case scenario")
            return
        end
        append!(solution[min_index], node)
    end
    return solution
end

function run_first_heuristic(inputFile::String, warmStart::Bool=false)::Union{Nothing, Dict{Int,Array}, Tuple{Float64, Float64}}
    start = time()
    println("Trying Mode 0")
    solution = construct_solution(inputFile)
    #println(solution)
    if !(check_feasibility(inputFile, solution))
        println("Trying Mode 1")
        solution = construct_solution(inputFile, 1)
        if !(check_feasibility(inputFile, solution))
            println("Trying Mode 2")
            solution = construct_solution(inputFile, 2)
            if !(check_feasibility(inputFile, solution))
                println("Trying Mode 3")
                solution = construct_solution(inputFile, 3)
                if !(check_feasibility(inputFile, solution))
                    println("Heuristic solution for ", inputFile, " is not feasible in the worst case scenario")
                    return
                end
            end
        end
    end
    if warmStart
        return solution
    else
        return time() - start, compute_worst_case(inputFile, solution)
    end
end

function run_heuristic(inputFile::String, timeLimit::Float64=-1.)::Union{Nothing,Tuple{Float64,Float64}}
    println("Solving instance ", inputFile, " with time limit " * string(timeLimit))
    start = time()
    best_solution = Dict{Int,Array}()
    best_res = -1
    best_mode = -1
    for mode in 0:2
        if timeLimit > 0 && (time() - start) > timeLimit
            println("No time left for different mode. Spent time is " * string(time() - start))
            break
        end
        solution = construct_solution(inputFile, mode)
        if (check_feasibility(inputFile, solution))
            res = compute_worst_case(inputFile, solution)
            if best_res == -1
                best_solution = solution
                best_res = res
                best_mode = mode
            elseif res < best_res
                best_solution = solution
                best_res = res
                best_mode = mode
            end
        end
    end
    if best_res == -1
        println("Heuristic solution for ", inputFile, " is not feasible in the worst case scenario")
        return
    end
    #println("Mode ", best_mode, " is the best one")
    leftTime = timeLimit - (time() - start)
    if timeLimit< 0
        solution = local_search_ter(inputFile, best_solution)
    else
        solution = local_search_ter(inputFile, best_solution, leftTime)
    end
    if !(check_feasibility(inputFile, solution))
        println("Local search result is not feasible for instance ", inputFile)
    end
    runTime = time() - start
    return runTime, compute_worst_case(inputFile, solution)
end

function run_all_instances()::Nothing
    for inputFile in DATA_FILES
        println("Solving instance ", inputFile)
        println(run_heuristic(inputFile))
    end
end

function find_seed(inputFile::String)::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    for k in 251:1000
        Random.seed!(k)
        solution = construct_solution(inputFile, 3)
        if (check_feasibility(inputFile, solution))
            println(k)
            return k
        end
    end
    println("rip")
    return
end

function local_search_one_iteration(inputFile::String, solution::Dict{Int,Array})::Bool
    include(DATA_DIR_PATH * "\\" * inputFile)
    best_score = compute_worst_case(inputFile, solution)
    best_switch = (0, 0, 0, 0)
    for (a, cluster_a) in solution
        for (b, cluster_b) in solution
            if a != b
                for i in 1:length(cluster_a)
                    for j in 1:length(cluster_b)
                        inter = cluster_a[i]
                        cluster_a[i] = cluster_b[j]
                        cluster_b[j] = inter
                        if (compute_worst_weight_sort(inputFile, cluster_a) < B) && (compute_worst_weight_sort(inputFile, cluster_b) < B)
                            new_score = compute_worst_case(inputFile, solution)
                            if (new_score < best_score)
                                best_score = new_score
                                best_switch = (a, b, i, j)
                            end
                        end
                        inter = cluster_a[i]
                        cluster_a[i] = cluster_b[j]
                        cluster_b[j] = inter
                    end
                end
            end
        end
    end
    if best_switch == (0, 0, 0, 0)
        return false
    else
        (a, b, i, j) = best_switch
        inter = solution[a][i]
        solution[a][i] = solution[b][j]
        solution[b][j] = inter
        return true
    end
end

function compute_static_distance_difference(distances::Any, cluster_a::Array, cluster_b::Array, i::Int, j::Int)::Float64
    remove = 0
    add = 0
    for k in 1:length(cluster_a)
        if k != i
            remove += distances[cluster_a[k], cluster_a[i]]
            add += distances[cluster_a[k], cluster_b[j]]
        end
    end
    for l in 1:length(cluster_b)
        if l != j
            remove += distances[cluster_b[l], cluster_b[j]]
            add += distances[cluster_b[l], cluster_a[i]]
        end
    end
    return add - remove
end

function compute_robust_distance_difference(inputFile::String, supplementary_distances::Vector{Any}, robust_cost::Int, cluster_a::Array, cluster_b::Array, i::Int, j::Int)::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    copy_supplementary_distances = copy(supplementary_distances)
    for k in 1:length(cluster_a)
        if k != i
            deleteat!(copy_supplementary_distances, indexin(lh[cluster_a[k]] + lh[cluster_a[i]], copy_supplementary_distances))
            append!(copy_supplementary_distances, lh[cluster_a[k]] + lh[cluster_b[j]])
        end
    end
    for l in 1:length(cluster_b)
        if l != j
            deleteat!(copy_supplementary_distances, indexin(lh[cluster_b[l]] + lh[cluster_b[j]], copy_supplementary_distances))
            append!(copy_supplementary_distances, lh[cluster_b[l]] + lh[cluster_a[i]])
        end
    end
    copy_supplementary_distances = reverse!(sort(copy_supplementary_distances))
    new_robust_cost = 0
    remaining_uncertainty = L
    for d in copy_supplementary_distances
        error = min(3, remaining_uncertainty)
        remaining_uncertainty -= error
        new_robust_cost += error * d
    end
    return new_robust_cost - robust_cost, copy_supplementary_distances
end

function local_search_one_iteration_bis(inputFile::String, distances::Any, solution::Dict{Int,Array}, supplementary_distances::Vector{Any}, robust_cost::Int, timeLimit::Float64=-1.)::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    start = time()
    best_static_gap = 0
    best_robust_gap = 0
    best_set_distances = []
    best_switch = (0, 0, 0, 0)
    for (a, cluster_a) in solution
        for (b, cluster_b) in solution
            if a != b
                for i in 1:length(cluster_a)
                    for j in 1:length(cluster_b)
                        if timeLimit > 0 && (time() - start) > timeLimit
                            return false, [], 0
                        end
                        inter = cluster_a[i]
                        cluster_a[i] = cluster_b[j]
                        cluster_b[j] = inter
                        if (compute_worst_weight_sort(inputFile, cluster_a) <= B) && (compute_worst_weight_sort(inputFile, cluster_b) <= B)
                            inter = cluster_a[i]
                            cluster_a[i] = cluster_b[j]
                            cluster_b[j] = inter
                            static_gap = compute_static_distance_difference(distances, cluster_a, cluster_b, i, j)
                            robust_gap, set_distances = compute_robust_distance_difference(inputFile, supplementary_distances, robust_cost, cluster_a, cluster_b, i, j)
                            if (static_gap + robust_gap < best_static_gap + best_robust_gap)
                                best_static_gap = static_gap
                                best_robust_gap = robust_gap
                                best_set_distances = set_distances
                                best_switch = (a, b, i, j)
                            end
                        else
                            inter = cluster_a[i]
                            cluster_a[i] = cluster_b[j]
                            cluster_b[j] = inter
                        end
                    end
                end
            end
        end
    end
    if best_switch == (0, 0, 0, 0)
        return false, [], 0
    else
        (a, b, i, j) = best_switch
        inter = solution[a][i]
        solution[a][i] = solution[b][j]
        solution[b][j] = inter
        robust_cost += best_robust_gap
        println(best_static_gap + best_robust_gap)
        return true, best_set_distances, robust_cost
    end
end

function local_search(inputFile::String, solution::Dict{Int,Array})::Dict{Int,Array}
    include(DATA_DIR_PATH * "\\" * inputFile)
    start = time()
    iteration_count = 0
    while local_search_one_iteration(inputFile, solution) && time() - start < 300
        iteration_count += 1
    end
    println("Number of Local_search iterations : ", iteration_count)
    return solution
end

function local_search_bis(inputFile::String, solution::Dict{Int,Array})::Dict{Int,Array}
    include(DATA_DIR_PATH * "\\" * inputFile)
    start = time()
    _, supplementary_distances, robust_cost, distances = compute_worst_case(inputFile, solution, 1)
    iteration_count = 1
    change = true
    while change && time() - start < 300
        iteration_count += 1
        change, supplementary_distances, robust_cost = local_search_one_iteration_bis(inputFile, distances, solution, supplementary_distances, robust_cost)
    end
    println("Number of Local_search iterations : ", iteration_count)
    return solution
end

function local_search_ter(inputFile::String, solution::Dict{Int,Array}, timeLimit::Float64=-1.)::Dict{Int,Array}
    include(DATA_DIR_PATH * "\\" * inputFile)
    println("Local search time limit is " * string(timeLimit))
    start = time()
    _, supplementary_distances, robust_cost, distances = compute_worst_case(inputFile, solution, 1)
    iteration_count = 1
    change = true
    runTime = time() - start
    while change && (timeLimit < 0 || runTime < timeLimit)
        start = time()
        iteration_count += 1
        change, supplementary_distances, robust_cost = local_search_one_iteration_bis(inputFile, distances, solution, supplementary_distances, robust_cost, timeLimit - runTime)
        runTime += time() - start
    end
    println("Number of Local_search iterations : ", iteration_count, " local search runtime " * string(runTime))
    return solution
end

function translate_output(inputFile::String, solution::Dict{Int,Array})::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    x = zeros(Int64, n, n)
    y = zeros(Int64, n, K)
    for (k, cluster_k) in solution
        for i in cluster_k
            y[i,k] = 1
            for j in cluster_k
                if i < j
                    x[i,j] = 1
                end
            end
        end
    end
    return x, y
end
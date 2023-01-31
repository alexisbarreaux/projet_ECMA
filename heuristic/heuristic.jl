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
            println("Worst weight : ", worst_weight, ", B : ", B)
            return false
        end
    end
    return true
end

function compute_worst_case(inputFile::String, solution::Dict{Int,Array})::Any
    include(DATA_DIR_PATH * "\\" * inputFile)
    l = computeDistances(coordinates)
    total_distance = 0
    supplementary_distances = []
    for (_, cluster) in solution
        for a in cluster
            for b in cluster
                if a != b
                    append!(supplementary_distances, lh[a] + lh[b])
                    total_distance += l[a, b]
                end
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

"""function run_heuristic(inputFile::String)::Any
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
    return compute_worst_case(inputFile, solution)
end"""

function run_heuristic(inputFile::String)::Any
    best_solution = Dict{Int,Array}()
    best_res = -1
    best_mode = -1
    for mode in 0:2
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
    println("Mode ", best_mode, " is the best one")
    return best_res
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
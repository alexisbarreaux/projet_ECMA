# This file contains code to build the results table and diagram
using DataFrames
using CSV

include("./utils/constants.jl")

"""
include("./results.jl")
resultsTable()
"""

function resultsTable()
    # Building base table
    numberOfInstances = length(DATA_FILES)
    results = DataFrame(instance=DATA_FILES, PR=[-1. for i in 1:numberOfInstances]
    , Dual_time=[-1. for i in 1:numberOfInstances], Dual_gap=["---" for i in 1:numberOfInstances]
    , Branch_time=[-1. for i in 1:numberOfInstances], Branch_gap=["---" for i in 1:numberOfInstances]
    , Plans_time=[-1. for i in 1:numberOfInstances], Plans_gap=["---" for i in 1:numberOfInstances]
    , Heur_time=[-1. for i in 1:numberOfInstances], Heur_gap=["---" for i in 1:numberOfInstances]
    , Heur_local_time=[-1. for i in 1:numberOfInstances], Heur_local_gap=["---" for i in 1:numberOfInstances]
    )

    # Loading needed files
    # Dual
    filePath =RESULTS_DIR_PATH * "\\" * "dual_5" * ".csv"
    dualRes = DataFrame(CSV.File(filePath))
    # B&cut
    filePath =RESULTS_DIR_PATH * "\\" * "branch_and_cut_5" * ".csv"
    branchRes = DataFrame(CSV.File(filePath))
    # Cuts
    filePath =RESULTS_DIR_PATH * "\\" * "cut_5" * ".csv"
    cutRes = DataFrame(CSV.File(filePath))
    # Heuristic first
    filePath =RESULTS_DIR_PATH * "\\" * "heuristic" * ".csv"
    heurRes = DataFrame(CSV.File(filePath))
    # Heuristic local search
    filePath =RESULTS_DIR_PATH * "\\" * "heuristic_local" * ".csv"
    heurLocalRes = DataFrame(CSV.File(filePath))

    for i in 1:nrow(results)
        instance = results[i, "instance"]
        
        # Get all rows
        dualRow = findfirst(==(instance), dualRes.instance)
        branchRow = findfirst(==(instance), branchRes.instance)
        cutRow = findfirst(==(instance), cutRes.instance)
        heurRow = findfirst(==(instance), heurRes.instance)
        heurLocalRow = findfirst(==(instance), heurLocalRes.instance)

        # Get best bound
        bestBound = -1.
        # Dual
        if dualRow != nothing
            if dualRes[dualRow, "gap"] != 1
                bestBound = dualRes[dualRow, "value"] / (dualRes[dualRow, "gap"] + 1)
            end
        end
        # B&cut
        if branchRow != nothing
            if branchRes[branchRow, "gap"] != 1
                branchBound = branchRes[branchRow, "value"] / (branchRes[branchRow, "gap"] + 1)
                bestBound = max(branchBound, bestBound)
            end
        end
        # Cuts
        if cutRow != nothing
            if cutRes[cutRow, "gap"] != 1
                cutBound = cutRes[cutRow, "value"] / (cutRes[cutRow, "gap"] + 1)
                bestBound = max(cutBound, bestBound)
            end
        end

        # Storing results
        # Dual
        if dualRow != nothing
            results[i, "Dual_time"] = dualRes[dualRow, "time"]
            if bestBound != -1.
                results[i, "Dual_gap"] = string(round(100*((dualRes[dualRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end

        # B&cut
        if branchRow != nothing
            results[i, "Branch_time"] = branchRes[branchRow, "time"]
            if bestBound != -1.
                results[i, "Branch_gap"] = string(round(100*((branchRes[branchRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end

        # Cuts
        if cutRow != nothing
            results[i, "Plans_time"] = cutRes[cutRow, "time"]
            if bestBound != -1.
                results[i, "Plans_gap"] = string(round(100*((cutRes[cutRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end

        # Heuristic first
        if heurRow != nothing
            results[i, "Heur_time"] = heurRes[heurRow, "time"]
            if bestBound != -1.
                results[i, "Heur_gap"] = string(round(100*((heurRes[heurRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end

        # Heuristic local
        if heurLocalRow != nothing
            results[i, "Heur_local_time"] = heurLocalRes[heurLocalRow, "time"]
            if bestBound != -1.
                results[i, "Heur_local_gap"] = string(round(100*((heurLocalRes[heurLocalRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end
    end
    # Saving to file
    resultsPath =RESULTS_DIR_PATH * "\\" * RESULTS_TABLE_FILE * ".csv"
    CSV.write(resultsPath, results, delim=";")

    return results
end
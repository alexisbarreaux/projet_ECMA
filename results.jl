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
    results = DataFrame(instance=DATA_FILES, PR=["---" for i in 1:numberOfInstances]
    , Dual_time=["---" for i in 1:numberOfInstances], Dual_gap=["---" for i in 1:numberOfInstances]
    , Branch_time=["---" for i in 1:numberOfInstances], Branch_gap=["---" for i in 1:numberOfInstances]
    , Plans_time=["---" for i in 1:numberOfInstances], Plans_gap=["---" for i in 1:numberOfInstances]
    , Heur_time=["---" for i in 1:numberOfInstances], Heur_gap=["---" for i in 1:numberOfInstances]
    , Heur_local_time=["---" for i in 1:numberOfInstances], Heur_local_gap=["---" for i in 1:numberOfInstances]
    )

    # Loading needed files
    # Static
    filePath =RESULTS_DIR_PATH * "\\" * "static_5" * ".csv"
    staticRes = DataFrame(CSV.File(filePath))
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
        staticRow = findfirst(==(instance), staticRes.instance)
        dualRow = findfirst(==(instance), dualRes.instance)
        branchRow = findfirst(==(instance), branchRes.instance)
        cutRow = findfirst(==(instance), cutRes.instance)
        heurRow = findfirst(==(instance), heurRes.instance)
        heurLocalRow = findfirst(==(instance), heurLocalRes.instance)

        # Get best bound and best feasible
        bestBound = -1.
        bestFeasible =1e63
        oneOpt = false
        # Dual
        if dualRow != nothing
            if dualRes[dualRow, "gap"] != 1
                bestBound = dualRes[dualRow, "value"] / (dualRes[dualRow, "gap"] + 1)
            end
            if !oneOpt && dualRes[dualRow, "optimal"]
                oneOpt = true
            end
        end
        # B&cut
        if branchRow != nothing
            if branchRes[branchRow, "gap"] != 1
                branchBound = branchRes[branchRow, "value"] / (branchRes[branchRow, "gap"] + 1)
                bestBound = max(branchBound, bestBound)
            end
            if !oneOpt && branchRes[branchRow, "optimal"]
                oneOpt = true
            end
        end
        # Cuts
        if cutRow != nothing
            if cutRes[cutRow, "gap"] != 1
                cutBound = cutRes[cutRow, "value"] / (cutRes[cutRow, "gap"] + 1)
                bestBound = max(cutBound, bestBound)
            end
            if !oneOpt && cutRes[cutRow, "optimal"]
                oneOpt = true
            end
        end
        
        
        # Storing results
        # Dual
        if dualRow != nothing
            results[i, "Dual_time"] = string(round(dualRes[dualRow, "time"], digits=2)) *"s"
            bestFeasible = dualRes[dualRow, "value"]
            if bestBound != -1.
                results[i, "Dual_gap"] = string(round(100*((dualRes[dualRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end
        
        # B&cut
        if branchRow != nothing
            results[i, "Branch_time"] = string(round(branchRes[branchRow, "time"], digits=2)) *"s"
            bestFeasible = min(branchRes[branchRow, "value"], bestFeasible)
            if bestBound != -1.
                results[i, "Branch_gap"] = string(round(100*((branchRes[branchRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end
        
        # Cuts
        if cutRow != nothing
            results[i, "Plans_time"] = string(round(cutRes[cutRow, "time"], digits=2)) *"s"
            bestFeasible = min(cutRes[cutRow, "value"], bestFeasible)
            if bestBound != -1.
                results[i, "Plans_gap"] = string(round(100*((cutRes[cutRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end

        # Heuristic first
        if heurRow != nothing
            results[i, "Heur_time"] = string(round(heurRes[heurRow, "time"], digits=2)) *"s"
            bestFeasible = min(heurRes[heurRow, "value"], bestFeasible)
            if bestBound != -1.
                results[i, "Heur_gap"] = string(round(100*((heurRes[heurRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end

        # Heuristic local
        if heurLocalRow != nothing
            results[i, "Heur_local_time"] = string(round(heurLocalRes[heurLocalRow, "time"], digits=2)) *"s"
            bestFeasible = min(heurLocalRes[heurLocalRow, "value"], bestFeasible)
            if bestBound != -1.
                results[i, "Heur_local_gap"] = string(round(100*((heurLocalRes[heurLocalRow, "value"]/bestBound) - 1), digits=1)) * "%"
            end
        end

        # PR
        # TODO faire une autre colonne quand on sait pas
        if staticRow != nothing
            prValue = string(round(100*(bestFeasible/staticRes[staticRow, "value"]  - 1), digits=1))
            if staticRes[staticRow, "optimal"]
                if oneOpt
                    results[i, "PR"] = prValue* "%"
                else    
                    results[i, "PR"] = "<= " * prValue * "%"
                end
            else
                results[i, "PR"] = "(?) " * prValue * "%"
            end
        end
    end
    # Saving to file
    resultsPath =RESULTS_DIR_PATH * "\\" * RESULTS_TABLE_FILE * ".csv"
    CSV.write(resultsPath, results, delim=";")

    return results
end
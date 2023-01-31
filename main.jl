using DataFrames
using CSV

include("./models/staticModel.jl")
include("./models/dualModel.jl")
include("./models/cutsModel.jl")
include("./models/branchAndCutModel.jl")

include("./utils/constants.jl")

"""
include("./main.jl")
staticSolveChosenInstance("30_eil_9.tsp", 5.0)
solveAllInstances()
solveAllInstances(cutSolve, 5.0, CUT_RESULTS_FILE)
solveAllInstances(brandAndCutSolve, 5.0, BRANCH_AND_CUT_RESULTS_FILE)
solveAllInstances(dualSolve, 30.0, DUAL_RESULTS_FILE)
staticSolveUnoptimizedInstance(60.0)
"""

function findUnoptimzedInstance(currentResults::DataFrame)::Tuple{String, Int}
    fileToRun = ""
    rowToReplace = -1
    for i in 1:nrow(currentResults)
        row = currentResults[i,:]
        if !row.optimal
            fileToRun = row.instance
            rowToReplace = i
            break
        end
    end
    return fileToRun, rowToReplace
end

function runInstanceAndUpdateDataframe(method::Function, currentResults::DataFrame, fileToRun::String, timeLimit::Float64, rowToReplace::Union{Int, Nothing}=nothing)::Bool
    
    result = method(fileToRun, false, true, timeLimit)
    if result == nothing
        println("NOT FEASIBLE!!")
        return false
    end

    optimal, solveTime, value, gap = result
    # Modify dataframe
    if rowToReplace == nothing
        rowToReplace = findfirst(==(fileToRun), currentResults.instance)
        if rowToReplace == nothing
            println("Pushing new row to results dataframe")
            push!(currentResults, [fileToRun optimal solveTime value gap])
            return true
        else
            currentRow= currentResults[rowToReplace,:]
            if value < currentRow.value || (solveTime > currentRow.time && value <= currentRow.value + 1e-6)
                println("Improved value for " * fileToRun)
                currentResults[rowToReplace,:] = [fileToRun optimal solveTime value gap]
                return true
            end
        end
    else
        currentRow = currentResults[rowToReplace,:]
        if value < currentRow.value || (solveTime > currentRow.time && value <= currentRow.value + 1e-6)
            println("Improved value for " * fileToRun)
            currentResults[rowToReplace,:] = [fileToRun optimal solveTime value gap]
            return true
        end
    end
    return false
end

function solveAllInstances(method::Function = staticSolve, timeLimit::Float64=-1., resultFile::String=STATIC_RESULTS_FILE)::Nothing
    # Loading
    filePath =RESULTS_DIR_PATH * "\\" * resultFile * ".csv"
    # Get unoptimal instance
    if !isfile(filePath)
        currentResults = DataFrame(instance=String[], optimal=Bool[], time=Float64[], value=Float64[], gap=Float64[])
    else
        currentResults = DataFrame(CSV.File(filePath))
    end

    # Run
    for fileToRun in ["318_lin_3.tsp", "318_lin_6.tsp", "318_lin_9.tsp",
        "400_rd_3.tsp", "400_rd_6.tsp", "400_rd_9.tsp",
        "532_att_3.tsp", "532_att_6.tsp", "532_att_9.tsp"]
        updatedDf = runInstanceAndUpdateDataframe(method, currentResults, fileToRun, timeLimit)
        if updatedDf
            CSV.write(filePath, currentResults, delim=";")
        end
    end
    return 
end

function solveAllNonOptimalInstances(method::Function = staticSolve, timeLimit::Float64=-1., resultFile::String=STATIC_RESULTS_FILE)::Nothing
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
    for fileToRun in DATA_FILES
        """
        if !endswith(fileToRun, "_3.tsp")
            continue
        end
        """
        rowIndex = findfirst(==(fileToRun), currentResults.instance)
        if rowIndex == nothing || !currentResults[rowIndex,:].optimal
            updatedDf = runInstanceAndUpdateDataframe(method, currentResults, fileToRun, timeLimit)
            if updatedDf
                CSV.write(filePath, currentResults, delim=";")
            end
            continue
        else
            println("Skipping " * fileToRun * ": already optimal.")
            continue
        end
    end
    return 
end

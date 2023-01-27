using DataFrames
using CSV

include("./models/staticModel.jl")
include("./utils/constants.jl")

"""
include("./main.jl")
staticSolveChosenInstance("30_eil_9.tsp", 5.0)
staticSolveAllInstances(60.0)
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

function staticRunInstanceAndUpdateDataframe(currentResults::DataFrame, fileToRun::String, timeLimit::Float64, rowToReplace::Union{Int, Nothing}=nothing)::Bool
    optimal, solveTime, value = staticSolve(fileToRun, false, true, timeLimit)
    # Modify dataframe
    if rowToReplace == nothing
        rowToReplace = findfirst(==(fileToRun), currentResults.instance)
        if rowToReplace == nothing
            println("Pushing new row to results dataframe")
            push!(currentResults, [fileToRun optimal solveTime value])
            return true
        else
            currentRow= currentResults[rowToReplace,:]
            if value < currentRow.value || (solveTime > currentRow.time && value <= currentRow.value + 1e-6)
                println("Improved value for " * fileToRun)
                currentResults[rowToReplace,:] = [fileToRun optimal solveTime value]
                return true
            end
        end
    else
        currentRow = currentResults[rowToReplace,:]
        if value < currentRow.value || (solveTime > currentRow.time && value <= currentRow.value + 1e-6)
            println("Improved value for " * fileToRun)
            currentResults[rowToReplace,:] = [fileToRun optimal solveTime value]
            return true
        end
    end
    return false
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
        if staticRunInstanceAndUpdateDataframe(currentResults, fileToRun, timeLimit, rowToReplace)
            CSV.write(filePath, currentResults)
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
    if staticRunInstanceAndUpdateDataframe(currentResults, fileToRun, timeLimit)
        CSV.write(filePath, currentResults)
    end

    return 
end

function staticSolveAllInstances(timeLimit::Float64=60.0, resultFile::String=STATIC_RESULTS_FILE)::Nothing
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
    updatedDf = false
    for fileToRun in DATA_FILES
        updateDf = updateDf || staticRunInstanceAndUpdateDataframe(currentResults, fileToRun, timeLimit)
    end
    if updatedDf
        CSV.write(filePath, currentResults)
    end
    return 
end

using DataFrames
using CSV

include("./models/staticModel.jl")
include("./models/dualModel.jl")
include("./models/cutsModel.jl")
include("./models/branchAndCutModel.jl")

include("./utils/constants.jl")

include("./heuristic/heuristic.jl")

"""
include("./main.jl")
staticSolveChosenInstance("30_eil_9.tsp", 5.0)
solveAllInstances()
solveAllInstances(staticSolve, 5.0, STATIC_RESULTS_FILE)
solveAllInstances(cutSolve, 5.0, CUT_RESULTS_FILE)
solveAllInstances(brandAndCutSolve, 5.0, BRANCH_AND_CUT_RESULTS_FILE)
solveAllInstances(dualSolve, 30.0, DUAL_RESULTS_FILE)
solveAllHeuristic()

staticSolveUnoptimizedInstance(60.0)

filePath =RESULTS_DIR_PATH * "\\" * "static" * ".csv"
staticRes = DataFrame(CSV.File(filePath))
filePath =RESULTS_DIR_PATH * "\\" * "static_5" * ".csv"
static5Res = DataFrame(CSV.File(filePath))
sort!(staticRes, [order(:instance)])
sort!(static5Res, [order(:instance)])
static5Res.value ./ staticRes.value
relativeDf = DataFrame(instance=staticRes.instance, optimal=staticRes.optimal, relative_difference=abs.(static5Res.value - staticRes.value) ./ staticRes.value)
filePath =RESULTS_DIR_PATH * "\\" * "static_relative_30_to_5" * ".csv"
CSV.write(filePath, relativeDf, delim=";")

filePath =RESULTS_DIR_PATH * "\\" * "dual" * ".csv"
dualRes = DataFrame(CSV.File(filePath))
filePath =RESULTS_DIR_PATH * "\\" * "dual_5" * ".csv"
dual5Res = DataFrame(CSV.File(filePath))
dualRelative = DataFrame(instance=String[], optimal=Bool[], relative_difference=Float64[])
for instance in dualRes.instance
    dualRow = findfirst(==(instance), dualRes.instance)
    dual5Row = findfirst(==(instance), dual5Res.instance)
    push!(dualRelative, [instance dualRes[dualRow, :].optimal abs(dualRes[dualRow, :].value - dual5Res[dual5Row, :].value)/dualRes[dualRow, :].value])
end
filePath =RESULTS_DIR_PATH * "\\" * "dual_relative_30_to_5" * ".csv"
CSV.write(filePath, dualRelative, delim=";")

filePath =RESULTS_DIR_PATH * "\\" * "dual_cplex_cuts" * ".csv"
dualCutRes = DataFrame(CSV.File(filePath))
filePath =RESULTS_DIR_PATH * "\\" * "dual_5" * ".csv"
dual5Res = DataFrame(CSV.File(filePath))
dualRelative = DataFrame(instance=String[], optimal=Bool[], difference=Float64[], gap_difference=Float64[])
for instance in dualCutRes.instance
    println(instance)
    dual5Row = findfirst(==(instance), dual5Res.instance)
    dualCutRow = findfirst(==(instance), dualCutRes.instance)
    if !(dual5Row==nothing)
        gapDiff= dual5Res[dual5Row, :].gap - dualCutRes[dualCutRow, :].gap
        push!(dualRelative, [instance dual5Res[dual5Row, :].optimal (dual5Res[dual5Row, :].value - dualCutRes[dualCutRow, :].value)/dual5Res[dual5Row, :].value gapDiff] )
    else
        push!(dualRelative, [instance dualCutRes[dualCutRow, :].optimal dualCutRes[dualCutRow, :].value -1.] )
    end
end

filePath =RESULTS_DIR_PATH * "\\" * "dual_cuts_comparison" * ".csv"
CSV.write(filePath, dualRelative, delim=";")

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
    for fileToRun in DATA_FILES
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


function solveAllHeuristic(resultFile::String=HEURISTICS_RESULTS_FILE)::Nothing
    # Loading
    filePath =RESULTS_DIR_PATH * "\\" * resultFile * ".csv"
    # Get unoptimal instance
    if !isfile(filePath)
        currentResults = DataFrame(instance=String[], time=Float64[], value=Float64[])
    else
        currentResults = DataFrame(CSV.File(filePath))
    end

    # Run
    for fileToRun in DATA_FILES
        updatedDf = runHeuristicAndUpdateDataframe(currentResults, fileToRun)
        if updatedDf
            CSV.write(filePath, currentResults, delim=";")
        end
    end
    return 
end

function runHeuristicAndUpdateDataframe(currentResults::DataFrame, fileToRun::String, rowToReplace::Union{Int, Nothing}=nothing)::Bool
    
    result = run_heuristic(fileToRun) # run_first_heuristic(fileToRun) 
    if result == nothing
        println("NOT FEASIBLE!!")
        return false
    end

    solveTime, value = result
    # Modify dataframe
    if rowToReplace == nothing
        rowToReplace = findfirst(==(fileToRun), currentResults.instance)
        if rowToReplace == nothing
            println("Pushing new row to results dataframe")
            push!(currentResults, [fileToRun solveTime value])
            return true
        else
            currentRow= currentResults[rowToReplace,:]
            if value < currentRow.value || (solveTime > currentRow.time && value <= currentRow.value + 1e-6)
                println("Improved value for " * fileToRun)
                currentResults[rowToReplace,:] = [fileToRun solveTime value]
                return true
            end
        end
    else
        currentRow = currentResults[rowToReplace,:]
        if value < currentRow.value || (solveTime > currentRow.time && value <= currentRow.value + 1e-6)
            println("Improved value for " * fileToRun)
            currentResults[rowToReplace,:] = [fileToRun solveTime value]
            return true
        end
    end
    return false
end



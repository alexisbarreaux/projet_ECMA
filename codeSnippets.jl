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

filePath =RESULTS_DIR_PATH * "\\" * "heuristic_best" * ".csv"
heurBestRes = DataFrame(CSV.File(filePath))
filePath =RESULTS_DIR_PATH * "\\" * "dual_5" * ".csv"
dual5Res = DataFrame(CSV.File(filePath))
dualHeurComparison = DataFrame(instance=String[], optimal=Bool[], difference=Float64[])
for instance in dual5Res.instance
    println(instance)
    dual5Row = findfirst(==(instance), dual5Res.instance)
    heurRow = findfirst(==(instance), heurBestRes.instance)
    if heurRow == nothing
        push!(dualHeurComparison, [instance dual5Res[dual5Row, :].optimal 1] )
    else
        push!(dualHeurComparison, [instance dual5Res[dual5Row, :].optimal (dual5Res[dual5Row, :].value - heurBestRes[heurRow, :].value)/dual5Res[dual5Row, :].value] )
    end
end

filePath =RESULTS_DIR_PATH * "\\" * "dual_heur_comparison" * ".csv"
CSV.write(filePath, dualHeurComparison, delim=";")

filePath =RESULTS_DIR_PATH * "\\" * "dual_warm" * ".csv"
dualWarmRes = DataFrame(CSV.File(filePath))
filePath =RESULTS_DIR_PATH * "\\" * "dual_5" * ".csv"
dual5Res = DataFrame(CSV.File(filePath))
dualRelative = DataFrame(instance=String[], optimal=Bool[], difference=Float64[], gap_difference=Float64[])
for instance in dual5Res.instance
    println(instance)
    dual5Row = findfirst(==(instance), dual5Res.instance)
    dualWarmRow = findfirst(==(instance), dualWarmRes.instance)
    if !(dualWarmRow==nothing)
        gapDiff= dual5Res[dual5Row, :].gap - dualWarmRes[dualWarmRow, :].gap
        push!(dualRelative, [instance dual5Res[dual5Row, :].optimal (dual5Res[dual5Row, :].value - dualWarmRes[dualWarmRow, :].value)/dual5Res[dual5Row, :].value gapDiff] )
    else
        push!(dualRelative, [instance dual5Res[dual5Row, :].optimal dual5Res[dual5Row, :].value -1.] )
    end
end

filePath =RESULTS_DIR_PATH * "\\" * "dual_warm_comparison" * ".csv"
CSV.write(filePath, dualRelative, delim=";")
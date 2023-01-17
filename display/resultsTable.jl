using PlotlyJS

include("../utils/constants.jl")
include("../utils/jsonUtils.jl")

function displayResultsTable()::Any
    staticResults = jsonLoadDictFromFile(RESULTS_DIR_PATH * "\\" *STATIC_RESULTS_FILE)

    return plot(table(
        header=attr(values=["A Scores", "B Scores"],
                    line_color="darkslategray",
                    fill_color="lightskyblue",
                    align="left"),
        cells=attr(
            values=[[100, 90, 80, 90], # 1st column
            [95, 85, 75, 95]], # 2nd column
            line_color="darkslategray",
            fill_color="lightcyan",
            align="left"
        )),Layout(width=500, height=500))
end
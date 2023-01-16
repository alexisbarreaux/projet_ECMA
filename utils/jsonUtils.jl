import JSON

function jsonDropToFile(filePath::String,data::Any)::Nothing
    jsonData = JSON.json(data)
    open(filePath, "w") do file
        write(file, jsonData)
    end
    return
end
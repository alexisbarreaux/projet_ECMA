function computeDistances(coordinates::Matrix{Float64})::Matrix{Float64}
    coordinatesSize = size(coordinates)
    n = coordinatesSize[1]
    distances = Matrix{Float64}(undef, n, n)
    for i in 1:n
        x_i, y_i= coordinates[i,:] 
        for j in i:n
            if i == j
                distances[i,j]=distances[j,i]= 0.
            else
                x_j, y_j = coordinates[j,:]
                currentDistance = sqrt( (x_i - x_j)^2 + (y_i - y_j)^2)
                distances[i,j]=distances[j,i]= currentDistance
            end
        end
    end

    return distances
end
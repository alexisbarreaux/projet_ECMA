using CPLEX

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

function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)
    # context_id == CPX_CALLBACKCONTEXT_CANDIDATE si le callback est
    # appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end
    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)
    # S’il n’y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end
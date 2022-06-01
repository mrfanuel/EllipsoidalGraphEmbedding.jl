@doc raw"""
    rename_com_unique(community)

rename the community labels by the smallest labels of integers.
# Arguments
- `community::Array{Int64,1}` array of integers

# Output
- `community::Array{Int64,1}` renamed array of integers

"""
function rename_com_unique(community::Array{Int64,1})

    N = length(community)
    community_reordered = vec(zeros(Int64, N, 1))

    n_c_best = length(unique(community))
    index = unique(community)
    for i = 1:N
        for l = 1:n_c_best
            if community[i] == index[l]
                community_reordered[i] = l
            end
        end
    end
    return community_reordered
end

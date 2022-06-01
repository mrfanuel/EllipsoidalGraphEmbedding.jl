

function random_generate(p::Array{Float64,1})

    r::Float64 = rand()
    N::Int64 = length(p)
    cum = zeros(size(p))

    cum[1] = p[1]

    elem::Int64 = 0

    for i = 2:N
        @inbounds cum[i] = cum[i-1] + p[i]
    end

    if r < cum[1]
        elem = 1
    else
        for i = 1:N-1
            if r > cum[i] && r < cum[i+1]
                elem = i + 1
            end
        end
    end

    return elem

end


function random_generate_vec(p::Array{Float64,1}, n_r::Int64)

    elem_vec::Array{Int64,1} = vec(zeros(n_r, 1))
    elem::Int64 = 0

    ok::Bool = false
    member::Bool = false
    elem_vec[1] = random_generate(p)

    for i = 2:n_r

        ok = false
        while ok == false
            elem = random_generate(p)
            member = false
            j = 1
            while j < i
                if elem_vec[j] == elem
                    member = true
                end
                j += 1
            end
            if member == false
                ok = true
            end
        end
        elem_vec[i] = elem

    end

    return elem_vec

end



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


# function MultiplyVecBlas(A::SparseMatrixCSC{Int64,Int64}, d::Array{Int64,2}, s::Float64, H::Array{Float64,2}, f::Float64)

#     N = size(H, 1)
#     r = size(H, 2)
#     X1 = zeros(N, r)
#     X1 = A * H - (d / s) * (d' * H) + f * H

#     return X1

# end


# function ReNameCom(community_n::Array{Int64,1})

#     community_n_n = zeros(Int64, size(community_n))
#     N_n = maximum(size(community_n))
#     k = 1
#     community_n_n[1] = 1
#     for i = 2:N_n
#         if community_n[i] != community_n[i-1]
#             k = k + 1

#         end
#         community_n_n[i] = k
#     end

#     return community_n_n

# end
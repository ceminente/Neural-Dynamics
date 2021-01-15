using ProgressBars
using DelimitedFiles
#l = 1
mu = 1
#e = 0.0001
N = 10^4
n_updates= 10^5

function compute_transition_probability(state, N, e, mu, l)
    #this function computes the unnormalized probability associated to each reaction
    inactive = (N - sum(state))
    P_eps = e*inactive
    P_p = (l*state*inactive)/N
    P_m = mu*state
    P = vcat( P_p, P_m, P_eps)
    return P
end

for l in 1:2, e in [0.1, 0.01, 0.001, 0.0001]
    # tau for first appereance of avalanche
    x = rand(1)[1]
    tau = (1/(e*N))*log(1/x)

    # initial conditions
    state = [1]
    t_size = [1]
    t_life = [tau]

    # final data on sizes (sizes) and duration (lifes)
    sizes = Int64[]
    lifes = Float64[]

    #history of actives
    active_story = [1]
    tau_story = [tau]



    for Gillespie_update in ProgressBar(1:n_updates)
        # unnormalized probability
        P = compute_transition_probability(state,N,e,mu,l)
        # normalization factor
        a0 = sum(P)
        # normalized probabilities
        P = cumsum(P)/a0
        #random  number
        r = rand(1)[1]
        #index of reaction
        idx = findfirst(x->x>r, P)
        #number of active avalanches
        k = length(state)
        #time before new reaction
        x = rand(1)[1]
        tau = (1/a0)*log(1/x)

        #UPDATES
        if (idx == 2k + 1)
            push!(state, 1)
            push!(t_size, 1)
            push!(t_life, 0)
        elseif (idx < k+1)
            state[idx] += 1
            t_size[idx] += 1
        else
            idx -= k
            state[idx] -=1
            if state[idx] == 0
                push!(sizes, t_size[idx])
                push!(lifes, t_life[idx])
                deleteat!(state, idx)
                deleteat!(t_life, idx)
                deleteat!(t_size, idx)
            end
        end
        # increase life of active states
        t_life = t_life .+ tau
        # save actives
        push!(active_story, sum(state))
        push!(tau_story, tau)
    end

    tau_story = cumsum(tau_story)

    npzwrite("active_story_$(e)_$(l).npy", active_story)
    npzwrite("tau_story_$(e)_$(l).npy", tau_story)
    npzwrite("final_time_$(e)_$(l).npy", lifes)
    npzwrite("final_size_$(e)_$(l).npy", sizes)
end

"""
lifes_r = map(x->round.(x), lifes)
u=sort(unique(lifes_r))
d= [count(x->x==i, lifes_r) for i in u]
scatter(u[2:end], [d[2:end], (u[2:end].^(-2)) .* 1000000] , xaxis=:log, yaxis=:log, label ="lifes", marker = (:hexagon, 3, 0.8, :green, stroke(0, 0.2, :black, :dot)))

a=sort(unique(sizes))
b= [count(x->x==i, sizes) for i in a]
scatter(a,[b, (a.^(-3/2)) .* 10000 ], xaxis=:log, yaxis=:log, label= "sizes", marker = (:hexagon, 3, 0.8, :purple, stroke(0, 0.2, :black, :dot)))

"""

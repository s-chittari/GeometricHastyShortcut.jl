using SparseArrays
using LinearAlgebra
using XLSX
using Arpack
using Plots
using CSV
using DelimitedFiles
using DataFrames
using Tables
using ProgressBars

###############################################################
####################### Gobal variables #######################
###############################################################

N = 8 # number of lattice sites
h = 1 # field strength
J = -0.5 # coupling strength
kb = 1
T = 1
C_barrier = 4
lambda_cutoff = 10
state_cutoff = 0.05

nn_mat = XLSX.readdata("/Users/suprajachittari/Documents/assisted assembly/d8_NN_forcode.xlsx","Sheet1","A1:C8") # nearest neighbor matrix (to be read in) 

###############################################################
########################## Functions ##########################
###############################################################

function clockwise_rot(occ::Array)
    clockwise_vec = [2, 3, 4, 1, 6, 7, 8, 5]
    occ_new = zeros(length(occ))

    for i=1:length((occ))
        occ_new[i] = clockwise_vec[occ[i]]
    end
    
    return Int.(occ_new)
end

function reflect_90_1(occ::Array)

    reflect_vec = [1, 5, 6, 2, 4, 8, 7, 3]

    occ_new = zeros(length(occ))

    for i=1:length(occ)
        occ_new[i] = reflect_vec[occ[i]]
    end

    return Int.(occ_new)
end

function reflect_90_2(occ::Array)

    reflect_vec = [1, 4, 8, 5, 2, 3, 7, 6]

    occ_new = zeros(length(occ))

    for i=1:length(occ)
        occ_new[i] = reflect_vec[occ[i]]
    end

    return Int.(occ_new)
end


function reflect_90_3(occ::Array)

    reflect_vec = [6, 2, 1, 5, 7, 3, 4, 8]

    occ_new = zeros(length(occ))

    for i=1:length(occ)
        occ_new[i] = reflect_vec[occ[i]]
    end

    return Int.(occ_new)
end

function reflect_90_4(occ::Array)

    reflect_vec = [3, 2, 6, 7, 4, 1, 5, 8]

    occ_new = zeros(length(occ))

    for i=1:length(occ)
        occ_new[i] = reflect_vec[occ[i]]
    end

    return Int.(occ_new)
end


function reflect_180(occ::Array)
    reflect = (mod.(occ .+ 4, 8))
    refl_new = replace(reflect, 0 => 8)
    return Int.(refl_new)
end

function symmetry_generation(v::Array)
    occ = findall(x->x==1, v) # indexes of occupied positions
    #print(occ)
    symmetries = zeros(N,(4*6))
    rotations = 3

    symmetries[:,1] = v
    indx = 2

    for i=1:rotations
        occ = findall(x->x==1, symmetries[:,indx-1])
        occ_new = clockwise_rot(occ)
        symmetries[occ_new, indx] .= 1
        #println(symmetries[:,indx])

        indx += 1
    end

    occ_new = reflect_90_1(findall(x->x==1, v))
    symmetries[occ_new,indx] .= 1
    indx += 1

    for i=1:(rotations)
        occ_new = clockwise_rot(findall(x->x==1, symmetries[:,indx-1]))
        symmetries[occ_new, indx] .= 1
        indx += 1
    end

    occ_new = reflect_90_2(findall(x->x==1, v))
    symmetries[occ_new,indx] .= 1
    indx += 1

    for i=1:rotations
        occ_new = clockwise_rot(findall(x->x==1, symmetries[:,indx-1]))
        symmetries[occ_new, indx] .= 1
        indx += 1
    end

    occ_new = reflect_90_3(findall(x->x==1, v))
    symmetries[occ_new,indx] .= 1
    indx += 1

    for i=1:rotations
        occ_new = clockwise_rot(findall(x->x==1, symmetries[:,indx-1]))
        symmetries[occ_new, indx] .= 1
        indx += 1
    end

    occ_new = reflect_90_4(findall(x->x==1, v))
    symmetries[occ_new,indx] .= 1
    indx += 1

    for i=1:rotations
        occ_new = clockwise_rot(findall(x->x==1, symmetries[:,indx-1]))
        symmetries[occ_new, indx] .= 1
        indx += 1
    end

    occ_new = reflect_180(findall(x->x==1, v))
    symmetries[occ_new,indx] .= 1
    indx += 1

    for i=1:rotations
        occ_new = clockwise_rot(findall(x->x==1, symmetries[:,indx-1]))
        symmetries[occ_new, indx] .= 1
        indx += 1
    end

    return symmetries
end

function E_diff(k::Int,v::Array) # defined for a 0 -> 1 transition
    # Identify any nearest neighbors
    nn_index = nn_mat[k,:] # indexes of nearest neighbors
    nn = v[nn_index] # spins of nearest neighbors
    n_f = sum(nn)
    E_diff = - h - sum(J * nn) # Compute energy difference

    return E_diff, Int(n_f)
end

function rate_ij(g::Array, j::Int, l::Int, E_diff::Float64, n_f::Int) # transition from state j to state l
    S_diff = log(g[j]) - log(g[l]) 
    return exp((E_diff/(2*kb*T))-(S_diff/2)+((C_barrier*n_f)/(kb*T)))
end

###############################################################
################# Transitions between states ##################
###############################################################

function state_transitions(c_bath)
    local rates, rates_columns, rates_rows

    g = Vector() # degeneracy
    c = Vector() # accepted configurations
    rates_rows = Vector()
    rates_columns = Vector()
    rates = Vector()
    occupancy = Vector()

    #############################################################
    ####################### Zero Occupancy ######################
    #############################################################

    push!(c, zeros(N))
    push!(g, 1)
    push!(occupancy, 0)

    ###############################################################
    ####################### Single Occupancy ######################
    ###############################################################

    c_single = zeros(N)
    c_single[1] = 1
    push!(c, c_single)
    push!(g, N)
    push!(occupancy, 1)

    # MATRIX LOOP

    for i=2:N # occupancy counter
        occs_minus = findall(x->x==(i-1), occupancy)
        
        for j in occs_minus
            # Find all positons = 0
            empty_pos = findall(x->x==0, c[j])

            for k in empty_pos
                # spin flip at kth position
                c_test = deepcopy(c[j])
                c_test[k] = 1

                # generate all configurations by symmetry rules - changing indexes
                symmetries_all = symmetry_generation(c_test)

                # remove duplicates in symmetry Vector
                symmetries_unique = unique(symmetries_all, dims=2)

                # search for each symmetry config in c vector
                traffic_light = true

                for l=1:length(c)
                    for i=1:length(symmetries_unique[1,:])
                        if symmetries_unique[:,i] == c[l]
                            # find transition rate
                            E_diff1, n_f = E_diff(k, symmetries_unique[:,i])

                            # assign matrix element - forward
                            push!(rates_rows, l) # ending state
                            push!(rates_columns, j) # initial state
                            push!(rates, c_bath*rate_ij(g, j, l, E_diff1, n_f))

                            # assign matrix element - reverse
                            push!(rates_rows, j) # ending state
                            push!(rates_columns, l) # initial state
                            push!(rates, rate_ij(g, l, j, -E_diff1, n_f))

                            traffic_light = false
                            break   
                        end
                    end

                    if traffic_light == false
                        break
                    end
                end

                if traffic_light == true
                    # pass new element to c (c_i -> new c)
                    push!(c, c_test)
                    push!(g, length(symmetries_unique[1,:]))
                    push!(occupancy, sum(c_test))

                    # find transition rate
                    E_diff2, n_f = E_diff(k, c_test)

                    # assign matrix element - forward
                    push!(rates_rows, length(c)) # final state
                    push!(rates_columns, j) # initial state
                    push!(rates, c_bath*rate_ij(g, j, length(c), E_diff2, n_f))

                    # assign matrix element - reverse
                    push!(rates_rows, j) # ending state
                    push!(rates_columns, length(c)) # initial state
                    push!(rates, rate_ij(g, length(c), j, -E_diff2, n_f))

                end
            end
        end
    end

    #############################################################
    ####################### Zero To One Transition ##############
    #############################################################

    E_diff0, n_f0 = E_diff(1, c[1])

    push!(rates_rows, 2)
    push!(rates_columns, 1)
    push!(rates, c_bath*rate_ij(g, 1, 2, E_diff0, n_f0))

    push!(rates_rows, 1)
    push!(rates_columns, 2)
    push!(rates, rate_ij(g, 2, 1, -E_diff0, n_f0))

    return rates_rows, rates_columns, rates
end

###############################################################
######### Matrix Construction and Eigendecomposition ##########
###############################################################
function mat_construction(rates_rows, rates_columns, rates)
    num_states = length(unique(rates_rows))
    rate_dense = zeros(num_states, num_states)

    for i=1:length(rates_rows)
        rate_dense[rates_rows[i],rates_columns[i]] = rates[i]
    end
    
    # Assign diagonal terms
    for i=1:length(rate_dense[1,:])
        rate_dense[i,i] = -sum(rate_dense[:,i])
    end

    return rate_dense
end

function eigendecomp(mat)
    evals, evecs = eigen(mat)

    if maximum(abs.(imag(evals))) > 10^(-8)
            println("Imaginary eigenvalue: ", maximum(abs.(imag(evals))))
    end
        
    evals = real(evals)
    evecs = real(evecs)
    
    evals = sort(evals, rev=true) # sort in descending order
    evecs = evecs[:, sortperm(evals)] # rearrange matching eval list
    
    norms = zeros(length(evals))
    
    #normalize the evecs
    for i in 1:length(evals)
        if evecs[1,i] == 0
            print("first term is zero","\n")
        end
        if evecs[1,i] > 0
            norms[i] = sum(abs.(evecs[:,i]))
        end
        if evecs[1,i] < 0
            norms[i] = -sum(abs.(evecs[:,i]))
        end
    
        evecs[:,i] /= norms[i]
    end

    sign_evec = 0

    # check sign 
    for i in 1:length(evals)
        for j in 1:length(evecs[:,i])
            if abs(evecs[j,i]) > 0.001
                sign_evec = sign(evecs[j,i])
                break
            end
        end
        evecs[:,i] *= sign_evec
    end

    return evals, evecs
end

###############################################################
########## Timescale separation, dimension reduction ##########
###############################################################
function timescale_separation(evals, evecs)
    chosen_eval = 0

    # choose eigenvalue set with appropriate timescale separation
    for i=2:(length(evals)-1)
        # println(abs(evals[i] / evals[i+1]))
        if abs(evals[i+1] / evals[i]) >= lambda_cutoff
            # println("The chosen eigenvalue set ends at eigenvalue number ", i)
            chosen_eval = i
            break
        end
    end

    # examine union of relevant supports
    support = Vector()
    if chosen_eval == 0
        println("ERROR: No spectral gap found.")
        return 0, 0
    else
        for i=2:chosen_eval
            push!(support, findall(x->x>=state_cutoff, abs.(evecs[:,i])))
        end
        support_list = vcat(support...) # concatenate arrays
        return evals[chosen_eval], sort(unique(support_list))
    end
end

function slow_modes_matrix(rate_matrix, eval_cutoff, support)
    # set timescale
    deltaT = -1/(eval_cutoff*1000)

    W = exp(rate_matrix*deltaT)

    for i=2:1000
        W = exp(rate_matrix*deltaT) * W
    end

    # print(W)

    # create empty matrix for slow modes
    n_s = length(support)
    slow_mat = zeros(n_s, n_s)

    n = length(rate_matrix[:,1])

    # construct slow modes matrix
    for i=1:n_s # starting states
        init = zeros(n)
        init[support[i]] = 1

        final_full = (W*init) ./ (deltaT*1000)

        for j=1:n_s
            slow_mat[j,i] = final_full[support[j]]
        end
    end

    # normalize the slow modes matrix
    for i=1:n_s
        slow_mat[i,i] = 0
        slow_mat[i,i] = -sum(slow_mat[:,i])
    end

    return slow_mat
end

function mapping_space(rate_matrix, eval_cutoff, support, evals, equilib)
    # initialize relevant timescales
    tau = -1/(eval_cutoff)
    dt = -1/(last(evals)*100000)
    # W = (I + rate_matrix*dt)^(tau/dt)
    W = exp(tau*rate_matrix)
    n = length(rate_matrix[:,1])
    n_s = length(support)

    # construct mapping matrix
    mapping = zeros(n_s, n)

    for i=1:n # starting state (full space)
        init = zeros(n)
        init[i] = 1
        
        final_full = W*init

        for j=1:n_s # ending state (reduced space)
            mapping[j,i] = final_full[support[j]]
        end
    end

    # normalize mapping matrix
    mapping_norm = zeros(n_s, n)
    for i=1:n
        mapping_norm[:,i] = (mapping[:,i] ./ sum(mapping[:,i]))*sum(equilib[support])
    end

    return W #mapping_norm
end

function recovery(p_vec, equilib, support)
    n = length(equilib)
    n_s = length(support)

    p_full = zeros(n)

    for i=1:n
        if i in support
            idx = findall(x->x==i, support)
            p_full[i] = p_vec[idx[1]]
        else
            p_full[i] = equilib[i]
        end
    end

    # normalization check
    if (sum(p_full) - 1) > 1E-3
        print("Error: recovered polymer is not normalized.")
    else
        return p_full
    end
end

###############################################################
################## Functions to evolve a state ################
###############################################################
function evolve_state(rate_mat, map, Q, p_init, time)
    p_exact = zeros(length(rate_mat[:,1]),length(time))
    p_approx = zeros(length(Q[:,1]),length(time))

    # mapping
    p_reduced = map*p_init

    for t in 1:length(time)
        # exact solution
        p_exact[:,t] = exp(rate_mat*time[t])*p_init

        # approximate solution
        p_approx[:,t] = exp(Q*time[t])*p_reduced
    end

    return p_exact, p_approx
end

###############################################################
################# Quantities of interest ######################
###############################################################
function KL_divergence(p,q)
    KL_div = 0

    for i = 1:length(p)
        KL_div += p[i]*log(p[i]/q[i])
    end

    return KL_div
end

function shannon_entropy(p)
    entropy = Vector()

    for j = 1:length(p[1,:])
        S = 0

        for i = 1:length(p[:,j])
            S += p[i,j]*log(p[i,j])
        end

        push!(entropy, -1*S)
    end

    return entropy
end

function rare_states(p, support, equilib1, equilib2)
    states_env1 = findall(x -> x >= state_cutoff, equilib1)
    states_env2 = findall(x -> x >= state_cutoff, equilib2)
    states = union(states_env1, states_env2)

    max_prob = 0
    max_state = 0

    for i = 1:length(p[1,:])
        prob, state = findmax(p[:,i])

        if (support[state] ∉ states)
            if prob > max_prob
                max_prob = prob
                max_state = support[state]
            end
        end
    end

    return max_state, max_prob
end

###############################################################
################## Functions for sanity check #################
###############################################################

function normalization_check(mat)
    for i in 1:length(mat[1,:])
        if abs(sum(mat[:,i])) > 1E-6
            println("Matrix is not normalized. Sum is ", sum(mat[i,:]))
            break
        end
    end
end


###############################################################
########################## Set generation #####################
###############################################################
c_vec_short = collect(0.1:0.5:3.6)
constrained_envs = length(c_vec_short)
c_vec = collect(0.1:0.5:7.1)
init_envs = length(c_vec)
second_envs = init_envs
set_size = constrained_envs*(init_envs^4)

mapping_matrices = Matrix{Float64}[]
# support_vector = Vector()
set = zeros(24, set_size)
histogram = zeros(set_size)
set_index = 1

env_pattern = zeros(set_size, 5+init_envs)

equilibs = zeros(24, init_envs)

for i = 1:init_envs
    # global set_index

    # initialize environment
    c_i = c_vec[i]

    # initialize matrix
    rates_rows, rates_columns, rates = state_transitions(c_i)
    rate_mat = mat_construction(rates_rows, rates_columns, rates)
    normalization_check(rate_mat)
    

    # eigendecomposition
    evals, evecs = eigendecomp(rate_mat)

    # timescale separation
    largest_eval, support = timescale_separation(evals, evecs)
    map = mapping_space(rate_mat, largest_eval, support, evals, evecs[:,1])
    push!(mapping_matrices, map)
    # push!(support_vector, support)

    equilibs[:,i] = evecs[:,1]

    if maximum(abs.(map*evecs[:,1] - evecs[:,1])) > 1E-10
        println("Error in mapping matrix")
    end
    # println(rate_mat*equilibs[:,i])
    # pass equilibrium distribution to the matrix
    # set[:,set_index] = evecs[:,1]
    # env_pattern[set_index, 1] = i
    # histogram[i] = 1
    # set_index += 1
end

# set epsilon here
epsilon = 1E-7 # dummy variable to be rechanged

# epsilon = epsilon 
print("Epsilon is ", epsilon)

# print()
# c1_then_2 = recovery(mapping_matrices[1]*set[:,2], set[:,1], support_vector[1])
# println(recovery(mapping_matrices[10]*c1_then_2, set[:,10], support_vector[10]))

# epsilon = 1
# println(epsilon)

pairs = Vector()
# for i = 1:constrained_envs
#     for j = 1:init_envs
#         for k = 1:init_envs
#             first = mapping_matrices[k]*mapping_matrices[j]*equilibs[:,i]
#             second = mapping_matrices[k]*mapping_matrices[j]*mapping_matrices[k]*mapping_matrices[j]*equilibs[:,i]
#             third = mapping_matrices[k]*mapping_matrices[j]*mapping_matrices[k]*mapping_matrices[j]*mapping_matrices[k]*mapping_matrices[j]*equilibs[:,i]

#             if maximum(abs.(second-third)) > 1E-2
#                 push!(pairs, [i,j,k])
#             end
#         end
#     end
# end
max = 1E-1
pair = Vector()
for i = 1:constrained_envs
    for j = 1:init_envs
        for k = 1:init_envs
            first = mapping_matrices[i]*mapping_matrices[k]*mapping_matrices[j]*equilibs[:,i]
            second = mapping_matrices[i]*mapping_matrices[j]*mapping_matrices[k]*equilibs[:,i]

            if maximum(abs.(second-first)) > max
                max = maximum(abs.(second-first))
                pair = [i,j,k]
            end
        end
    end
end
print(pair)
print(max)
println(equilibs[:,1])
println(mapping_matrices[15]*mapping_matrices[3]*equilibs[:,1])
println(mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*equilibs[:,1])
println(mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*equilibs[:,1])
println(mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*equilibs[:,1])

println(mapping_matrices[3]*mapping_matrices[6]*mapping_matrices[15]*mapping_matrices[3]*equilibs[:,1])

println(equilibs[:,1])
println(mapping_matrices[3]*equilibs[:,1])
println(mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*equilibs[:,1])
println(mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*equilibs[:,1])
println(mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*mapping_matrices[15]*mapping_matrices[3]*equilibs[:,1])

println(mapping_matrices[8]*mapping_matrices[1]*mapping_matrices[15]*equilibs[:,8])
println(mapping_matrices[8]*mapping_matrices[15]*mapping_matrices[1]*equilibs[:,8])

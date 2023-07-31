const N = 8 # number of lattice sites
const h = 1 # field strength
const J = -0.5 # coupling strength
const kb = 1
const T = 1
const C_barrier = 4
const lambda_cutoff = 10
const state_cutoff = 0.05

const nn_mat = [4 2 5
                1 3 6
                2 4 7
                3 1 8
                8 6 1
                5 7 2
                6 8 3
                7 5 4]

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

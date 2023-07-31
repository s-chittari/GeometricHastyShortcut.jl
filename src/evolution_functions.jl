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

function normalization_check(mat)
    for i in 1:length(mat[1,:])
        if abs(sum(mat[:,i])) > 1E-6
            println("Matrix is not normalized. Sum is ", sum(mat[i,:]))
            break
        end
    end
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
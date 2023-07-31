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

    equilibs[:,i] = evecs[:,1]

    if maximum(abs.(map*evecs[:,1] - evecs[:,1])) > 1E-10
        println("Error in mapping matrix")
    end
end

# set epsilon here
epsilon = 1E-7 # dummy variable to be rechanged

# epsilon = epsilon 
print("Epsilon is ", epsilon)

for i = 1:init_envs
    for j = 1:init_envs
        for k = 1:init_envs
            global set_index

            dist = mapping_matrices[k]*mapping_matrices[j]*equilibs[:,i]

            set[:,set_index] = dist
            env_pattern[set_index, 1] = i
            env_pattern[set_index, 2] = j 
            env_pattern[set_index, 3] = k 

            for m = 1:init_envs
                mapped_dist = mapping_matrices[m]*dist
        
                if maximum(abs.(mapped_dist - equilibs[:,m])) < epsilon
                    env_pattern[set_index,(m+3)] = m
                end
            end

            set_index += 1
        end
    end
end

# generation of intermediate set

intermediate_set = zeros(24, init_envs^2)
intermediate_env = zeros(init_envs^2, 2)
new_index = 1

for i = 1:init_envs
    for j = 1:init_envs
        global new_index

        intermediate_set[:,new_index] = mapping_matrices[j]*equilibs[:,i]
        intermediate_env[new_index, 1] = i
        intermediate_env[new_index, 2] = j 

        new_index += 1
    end
end

equilib_df = DataFrame(equilibs, :auto)
CSV.write("equilib.csv", equilib_df)

set_df = DataFrame(set, :auto)
CSV.write("set.csv", set_df)

int_env_df = DataFrame(intermediate_env, :auto)
CSV.write("intermediate_env.csv", int_env_df)

int_set = DataFrame(intermediate_set, :auto)
CSV.write("intermediate_set.csv", int_set)

env_df = DataFrame(env_pattern, :auto)
CSV.write("environment_pattern.csv", env_df)
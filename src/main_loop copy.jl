
function main_loop_kernel(
    solutions,                        # [n_LPs, n_vars−1] = Holds solutions to LPs
    objectives,                       # [n_LPs] = Holds dual objectives
    original_variable_lower_bounds,   # [n_LPs, n_vars] = Used for getting unscaled convergence info
    original_variable_upper_bounds,   # [n_LPs, n_vars] = ""
    original_right_hand_side,         # [n_LPs × total_LP_length] = RHS h before scaling
    original_objective_vector,        # [n_LPs, n_vars] = Cost vector c before scaling of each LP objective
    original_objective_constant,      # [n_LPs] = the offset constant in the objective function of each LP
    # =============================================================================================================================#
    # The original_* versions are only used to compute unscaled convergence statistics
    # the algorithm runs entirely on the scaled_* versions.
    # =============================================================================================================================#
    scaled_variable_lower_bounds,     # [n_LPs, n_vars] = Lower bounds after initial rescaling using Ruiz/Pock-Chambolle rescaling
    scaled_variable_upper_bounds,     # [n_LPs, n_vars] = Upper bounds after initial rescaling using Ruiz/Pock-Chambolle rescaling
    scaled_constraint_matrix,         # [n_LPS × total_LP_length, n_vars] = Constraint matrix G~ after initial rescaling
    scaled_right_hand_side,           # [n_LPs × total_LP_length] = Right-hand h~ side after initial rescaling
    scaled_objective_vector,          # [n_LPs, n_vars] = Objective vector c~ after initial rescaling
    nz_count,                         # Int = Total number of nonzeros in the constraint matrix
    nz_rows,                          # [nz_count] = Row indices containing nonzeros in the constraint matrix
    nz_cols,                          # [nz_count] = Column indices containing nonzeros in the constraint matrix
    active_constraint,                # [n_LPs × total_LP_length] = Flag to indicate whether a constraint is in use.
    variable_rescaling,               # [n_LPs, n_vars] = Total value used for variable rescaling - Per-variable scale factor Dx so that unscaled_x = scaled_x / D_x
    constraint_rescaling,             # [n_LPs × total_LP_length] = Per-constraint scale factor Dc - ​Total value used for constraint rescaling
    # input_current_primal_solution,    # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_dual_solution,      # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_primal_product,     # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_dual_product,       # Usually 0s, but used for hot-starts (if desired in the future)
    # Last Restart Solutions
    # Current Iterate State - These are updated every PDLP step and represent where the algorithm currently is
    current_primal_solution,          # [n_LPs, n_vars] = Current x~ (scaled) - Currently active solution state
    current_dual_solution,            # [n_LPs × total_LP_length] = Current y~ (scaled)
    current_dual_product,             # [n_LPs, n_vars] = Cached G~^T * y~ to avoid recomputing
    current_primal_product,           # [n_LPs × total_LP_length] = Cached G~ * x~ to avoid recomputing (scaled)
    current_primal_gradient,           # [n_LPs, n_vars] = 
    initial_primal_solution,
    initial_dual_solution,
    pdhg_primal_solution,
    pdhg_dual_solution,
    # Unscaled Solution for convergence checking
    original_primal_solution,         # [n_LPs, n_vars] = Unscaled solution state information, used to calculate actual infeasibility
    original_primal_gradient,         # [n_LPs, n_vars] = Unscaled primal gradient
    original_dual_solution,           # [n_LPs × total_LP_length] = Unscaled dual solution y
    original_primal_product,          # [n_LPs × total_LP_length] = G~ * x~ unscaled (== G * x)
    # Step Computation Buffers (scratch space per PDLP step)
    buffer_kkt_primal_solution,       # [n_LPs, n_vars] = normalized x~ with primal_ray_norm for infeasibility ray checks Intermediate storage for calculating infeasibility
    buffer_kkt_primal_product,        # [n_LPs × total_LP_length] = normalized original_primal_product using primal_ray_norm (G~ * buffer_kkt_primal_solution)
    buffer_kkt_lower_variable_violation, # [n_LPs, n_vars] = max(l - x, 0) : Infeasible variable/constraint violation storage
    buffer_kkt_upper_variable_violation, # [n_LPs, n_vars] = max(x - u, 0) : Infeasible variable/constraint violation storage
    buffer_kkt_reduced_costs,         # [n_LPs, n_vars] = 
    delta_primal,                     # [n_LPs, n_vars] = Δx~ = x'~ - x~ : the value used in calculating PDLP steps
    delta_primal_product,             # [n_LPs × total_LP_length] = G~ * Δx~
    delta_dual,                       # [n_LPs × total_LP_length] = Δy~ = y'~ - y~ : the value used in calculating PDLP steps
    # Parameters for Control & Termination
    input_primal_weight,              # [n_LPs] = The original primal weight (could be moved internally if needed)
    input_step_size,                  # [n_LPs] = The original step size η (could be moved internally if needed)
    termination_reason,               # [n_LPs] = TerminationReason enum per LP : Field to let the user know why the LP terminated
    current_LP_length,                # Int = The current length of each individual LP (i.e., number of constraints)
    total_LP_length,                  # Int = The total possible size of each individual LP (i.e., max allowed number of constraints)
    n_LPs,                            # Int = The current number of LPs being solved (fewer than the max may be used during B&B)
    n_vars,                           # Int = The number of primal variables in the problem
    iteration_limit,                  # Int = The maximum number of PDLP steps allowed before termination
    kkt_matrix_pass_limit,            # Float = The maximum number of KKT matrix passes (generally unused)
    termination_evaluation_frequency, # Number of PDLP steps to take before checking termination criteria (default: 200)
    necessary_reduction_for_restart,  # Float = β_necessary = Necessary KKT error decrease factor needed for a restart if no progress is being made (default: 0.5)
    sufficient_reduction_for_restart, # Float = β_sufficient = Sufficient KKT error decrease factor that would trigger a restart immediately (default: 0.2)
    artificial_ratio_for_restart,     # Float = β_artificial = artificial restart to avoid long inner loop (default = )
    extrapolation_coefficient,        # Float = Value used in PDLP steps (default: 1.0)
    reflection_coefficient,           # Float = Reflection coefficient γ ∈ [0,1] at reflection step (default: 1.0)
    pid_KP,                           # Float = Primal Weight Update PID controller proportional coefficient (default: 0.99)
    pid_KI,                           # Float = Primal Weight Update PID controller integral coefficient (default: 0.01)
    pid_KD,                           # Float = Primal Weight Update PID controller derivative coefficient (default: 0)
    abs_tol,                          # Float = Absolute tolerance for termination
    rel_tol,                          # Float = Relative tolerance for termination
    eps_primal_infeasible,            # Float = Primal infeasibility tolerance
    eps_dual_infeasible,              # Float = Dual infeasibility tolerance
    return_code,                      # Int = Indicator for returning primal obj (1), dual obj (2), or both (3)
    global_upper_bound,               # Float = Information about the B&B upper bound (PDLP terminates if a dual feasible solution is above this value)
    skip_hard_problems,               # Bool = Flag to skip problems with too many iterations
    global_counter,                   # [1] (atomic Int32) = Count of LPs being solved - Running count of successfully completed LPs (shared across blocks via atomics)
    iteration_counter,                # [1] (atomic Int32) = Count of total iterations for solved LPs - Total iterations over all completed LPs (used to compute average for skip_hard_problems)
    skip_flag,                        # [n_LPs] = Flag to completely skip an individual LP - Per-LP flag to skip entirely (e.g., already pruned in B&B)
    iterations,                       # [n_LPs] = The final number of iterations needed for each LP
    )
    # In this kernel, assume that the number of blocks is equal to the number
    # of LPs, so that all threads are working on one LP
    LP = blockIdx().x
    idx = threadIdx().x
    
    block_stride = blockDim().x # total number of threads in the block
    grid_stride = gridDim().x # total number of blocks in the grid

    # Calculate strides for parallel reductions as the largest power of 2
    # less than n_vars (or current_LP_length)
    var_stride = Int32(1) << floor(Int32, log2(n_vars))
    len_stride = Int32(1) << floor(Int32, log2(current_LP_length))

    # Set up dynamic shared space
    shared_space = @cuDynamicSharedMem(Float64, max(n_vars, current_LP_length))

    while LP <= n_LPs
        # Check if we're supposed to skip this LP
        if skip_flag[LP]
            LP += grid_stride
            continue
        end

        # Initialize basic information for this LP
            ## information bookkeeping
            ## Only thread 1 of this block will store these information 
            ## the information below does not need to be stored per thread
        if idx==1
            cumulative_kkt_passes = 0.5

            # Initialize temporary Float64 values needed for calculation
            last_restart_primal_distance_moved = 0.0
            last_restart_dual_distance_moved = 0.0
            restart_primal_distance = 0.0
            restart_dual_distance = 0.0
            buffer_kkt_dual_objective = 0.0
            buffer_kkt_dual_res_inf = 0.0
            CI_primal_objective = 0.0 # Needed for optimality termination check
            CI_dual_objective = 0.0
            CI_l2_primal_residual = 0.0 # Needed for optimality termination check
            CI_l2_dual_residual = 0.0
            last_restart_length = 1.0
            last_reduction_ratio = 1.0
            restart_error = 0.0
            sum_restart_error = 0.0
            last_restart_error = 0.0
            current_kkt_residual = 0.0
            candidate_kkt_residual = 0.0
            restart_choice = RESTART_CHOICE_NO_RESTART
            cross_term = 0
            squared_delta_dual = 0
            squared_delta_primal = 0
            anchor_cross_term = 0
            anchor_squared_delta_primal = 0
            anchor_squared_delta_dual = 0
            initial_fixed_error = 0.0
            candidate_fixed_error = 0.0
            old_error = 0.0
        end

        # Other information that every thread needs (stored in the Registery memory)
            ## The values are always identical across threads — it is replicated rather than shared for performance purposes.
        
        total_iterations = Int32(0) # total iteration counter (T)
        inner_iterations = Int32(0)
        

        primal_ray_norm = 0.0


        
        # Information that is much easier to save as static shared memory (L1 Block memory)
        do_restart = @cuStaticSharedMem(Bool, 1)
        primal_weight = @cuStaticSharedMem(Float64, 1)
        numerical_error = @cuStaticSharedMem(Bool, 1)
        step_size = @cuStaticSharedMem(Float64, 1)
        if idx==1
            do_restart[1] = false
            primal_weight[1] = input_primal_weight[LP]
            step_size[1] = input_step_size[LP]
            # step_size[1] = 0.998 #TODO: remove this later
            numerical_error[1] = false
        end
        
        # Set up the starting row for this LP (minus 1, so that the first
        # row to consider is `active_row + 1`)
        active_row = (LP-Int32(1)) * total_LP_length 

        # Set up current values to match inputs
        while idx <= n_vars
            current_primal_gradient[LP, idx] = scaled_objective_vector[LP, idx]
            idx += block_stride
        end
        # resetting idx to its original thread id after tweaking it to work on different variables above
        idx = threadIdx().x
        

        # Set up objective vector and RHS norms
        while idx <= n_vars
            shared_space[idx] = original_objective_vector[LP, idx]^2
            idx += block_stride
        end
        idx = threadIdx().x
        parallel_sum(shared_space, block_stride, var_stride, n_vars)

        if idx==1
            cache_l2_norm_primal_linear_objective = sqrt(shared_space[1])
        end

        while idx <= current_LP_length
            shared_space[idx] = original_right_hand_side[active_row + idx]^2
            idx += block_stride
        end
        idx = threadIdx().x
        parallel_sum(shared_space, block_stride, len_stride, current_LP_length)

        if idx==1
            cache_l2_norm_primal_right_hand_side = sqrt(shared_space[1])
        end

        # Begin the main loop
        ITERATION_COUNT = 0
        while total_iterations <= iteration_limit

                ##########################################################################################
                #                       Taking Halpern Reflected PDHG steps                              #
                ##########################################################################################

                epoch_iterations = Int32(0) # counter for inner loop for computation purposes

                while epoch_iterations < termination_evaluation_frequency

                    ITERATION_COUNT += 1
                    idx==1 && CUDA.@cuprintln(">>>>>>>>>>>>ITERATION $ITERATION_COUNT HAS BEGUN<<<<<<<<<<<<<<<<")

                    idx==1 && CUDA.@cuprintln("\nINITIAL VALUES")
                    idx==1 && CUDA.@cuprintln("Current primal solution: [$(current_primal_solution[1,1]), $(current_primal_solution[1,2])]")
                    idx==1 && CUDA.@cuprintln("Current primal product: [$(current_primal_product[1]), $(current_primal_product[2]), $(current_primal_product[3]), $(current_primal_product[4]), $(current_primal_product[5]), $(current_primal_product[6])]")
                    idx==1 && CUDA.@cuprintln("current_dual_product: [$(current_dual_product[1,1]), $(current_dual_product[1,2])]")
                    idx==1 && CUDA.@cuprintln("Current dual solution: [$(current_dual_solution[1]), $(current_dual_solution[2]), $(current_dual_solution[3]), $(current_dual_solution[4]), $(current_dual_solution[5]), $(current_dual_solution[6])]\n")
                    
                    local_primal_weight = unsafe_load(CUDA.pointer(primal_weight, 1))
                    local_step_size = unsafe_load(CUDA.pointer(step_size, 1))
                    
                    # Step 1) Add one to the total iterations tracker and cumulative kkt passes
                    
                    if idx==1
                        cumulative_kkt_passes += Int32(1)
                    end
                        
                    # Step 2) Operate on matrices to:
                    # A) Update pdhg_primal_solution (next_x= T(x[k])) -> pdhg_primal_solution will be used to initialze the anchor if restart is triggered

                    while idx <= n_vars 
                        pdhg_primal_solution[LP, idx] = min(scaled_variable_upper_bounds[LP, idx], max(scaled_variable_lower_bounds[LP, idx], current_primal_solution[LP, idx] - 
                                                        (local_step_size/local_primal_weight) * (scaled_objective_vector[LP, idx] - current_dual_product[LP, idx])))
                        idx += block_stride
                    end

                    idx = threadIdx().x

                    # Compute delta primal for dual step
                    while idx <= n_vars 
                        delta_primal[LP, idx] = pdhg_primal_solution[LP, idx] - current_primal_solution[LP, idx]
                        idx += block_stride
                    end

                    idx = threadIdx().x

                    # Reset delta_primal_product and the shared space
                    while idx <= current_LP_length
                        delta_primal_product[active_row + idx] = 0.0
                        shared_space[idx] = 0.0
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    sync_threads()

                    # Loop over nonzeros to update delta_primal_product
                    while idx <= nz_count
                        if active_constraint[active_row + nz_rows[idx]]
                            CUDA.atomic_add!(CUDA.pointer(shared_space, nz_rows[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * delta_primal[LP, nz_cols[idx]])
                        end
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    sync_threads()
                    while idx <= current_LP_length
                        delta_primal_product[active_row + idx] += shared_space[idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # Compute other terms normally
                    while idx <= current_LP_length
                        pdhg_dual_solution[active_row + idx] = max(0.0, (current_dual_solution[active_row + idx] + 
                                            (local_primal_weight * local_step_size) * (scaled_right_hand_side[active_row + idx] -
                                            ((Int32(1) + extrapolation_coefficient) * 
                                            delta_primal_product[active_row + idx]) - 
                                            (extrapolation_coefficient * current_primal_product[active_row + idx]))))

                        idx += block_stride
                    end
                    
                    idx = threadIdx().x

                    
                    # taking the Halpern scheme with Reflection step and calculate it as delta_primal
                    while idx <= n_vars
                        delta_primal[LP, idx] = ((inner_iterations + 1)/ (inner_iterations + 2)) * ( ( 1 + reflection_coefficient) * pdhg_primal_solution[LP, idx] - reflection_coefficient * current_primal_solution[LP, idx]) + (1 / (inner_iterations + 2)) * initial_primal_solution[LP, idx] - current_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    

                    while idx <= n_vars
                        current_primal_solution[LP, idx] += delta_primal[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    if idx == 1 
                        CUDA.@cuprintln("step = $(epoch_iterations): current primal= $(current_primal_solution[1,1]), $(current_primal_solution[1,2])")
                        # CUDA.@cuprintln("----------: $(pdhg_primal_solution[1,1]), $(pdhg_primal_solution[1,2])")
                    end

                    while idx <= current_LP_length

                        delta_dual[active_row + idx] = ((inner_iterations + 1)/ (inner_iterations + 2)) * ( ( 1 + reflection_coefficient) * pdhg_dual_solution[active_row + idx] - reflection_coefficient * current_dual_solution[active_row + idx]) + (1 / (inner_iterations + 2)) * initial_dual_solution[active_row + idx] - current_dual_solution[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    while idx <= current_LP_length
                        current_dual_solution[active_row + idx] += delta_dual[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    sync_threads()

                    # Reset delta_primal_product and the shared space
                    while idx <= current_LP_length
                        delta_primal_product[active_row + idx] = 0.0
                        shared_space[idx] = 0.0
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    sync_threads()

                    # Loop over nonzeros to update delta_primal_product
                    while idx <= nz_count
                        if active_constraint[active_row + nz_rows[idx]]
                            CUDA.atomic_add!(CUDA.pointer(shared_space, nz_rows[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * delta_primal[LP, nz_cols[idx]])
                        end
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    sync_threads()
                    while idx <= current_LP_length
                        delta_primal_product[active_row + idx] += shared_space[idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    
                    while idx <= current_LP_length
                        current_primal_product[active_row + idx] += delta_primal_product[active_row + idx]
                        shared_space[idx] = delta_primal_product[active_row + idx] * delta_dual[active_row + idx]
                        idx += block_stride
                    end

                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                   
                    if idx==1
                        cross_term = shared_space[1]
                        if inner_iterations == Int32(0)
                            anchor_cross_term  = cross_term
                        end
                    end

                    # Compute the l2 norm of delta_primal and delta_dual for CURRENT solution and if we are the at first inner iteration, storing them for anchor
                    while idx <= n_vars 
                        shared_space[idx] = delta_primal[LP, idx] ^ 2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx == 1
                        squared_delta_primal = shared_space[1]
                        if inner_iterations == Int32(0)
                            anchor_squared_delta_primal = squared_delta_primal
                        end
                    end

                    while idx <= current_LP_length
                        shared_space[idx] = delta_dual[active_row + idx] ^ 2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)

                    if idx == 1
                        squared_delta_dual = shared_space[1]
                        if inner_iterations == Int32(0)
                            anchor_squared_delta_dual = squared_delta_dual
                        end
                    end

                    # Reset current_dual_product and the shared space
                    while idx <= n_vars
                        current_dual_product[LP, idx] = 0.0
                        shared_space[idx] = 0.0
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    sync_threads()

                    # Loop over nonzeros to update current_dual_product: 
                    while idx <= nz_count
                        if active_constraint[active_row + nz_rows[idx]]
                            CUDA.atomic_add!(CUDA.pointer(shared_space, nz_cols[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * current_dual_solution[active_row + nz_rows[idx]])
                        end
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    sync_threads()
                    while idx <= n_vars
                        current_dual_product[LP, idx] += shared_space[idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    sync_threads() 


                    if inner_iterations == Int32(0)
                        if unsafe_load(CUDA.pointer(do_restart, 1)) == true 
                            if idx == 1
                                initial_fixed_error = sqrt((local_primal_weight / local_step_size) * anchor_squared_delta_primal + (1 / (local_step_size * local_primal_weight)) * anchor_squared_delta_dual + 2 * local_step_size * anchor_cross_term)
                                # initial_fixed_error = sqrt((local_primal_weight) * anchor_squared_delta_primal + (1 / (local_primal_weight)) * anchor_squared_delta_dual + 2 * local_step_size * anchor_cross_term)
                                CUDA.@cuprintln("initial fixed point error = $initial_fixed_error \n")
                                do_restart[1] = false
                            end
                        end
                    end

                    inner_iterations += Int32(1)
                    epoch_iterations += Int32(1)           
                end         

                sync_threads()
                # Update Current Primal Gradient 

                while idx <= n_vars
                    current_primal_gradient[LP, idx] = scaled_objective_vector[LP, idx] - current_dual_product[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x

                # # printing the progress 
                # if idx == 1
                #     CUDA.@cuprintln(": primal solution = [$(current_primal_solution[1,1]), $(current_primal_solution[1,2])]")
                # end
                ##########################################################################################
                #                            Compute Fixed-point error                                   #
                ##########################################################################################
                if idx==1

                    
                    # initial_fixed_error = sqrt((primal_weight[1]) * anchor_squared_delta_primal + (1 / (primal_weight[1])) * anchor_squared_delta_dual + 2 * anchor_cross_term)
                    candidate_fixed_error = sqrt((primal_weight[1] / step_size[1]) * squared_delta_primal + (1 / (step_size[1] * primal_weight[1])) * squared_delta_dual + step_size[1] * 2 * cross_term)
                    # candidate_fixed_error = sqrt((primal_weight[1]) * squared_delta_primal + (1 / (primal_weight[1])) * squared_delta_dual + step_size[1] * 2 * cross_term)
                    # CUDA.@cuprint("fixed point error = $(candidate_fixed_error)")
                    # candidate_fixed_error = sqrt((primal_weigh[1]) * squared_delta_primal + (1 / (primal_weight[1])) * squared_delta_dual + 2 * cross_term)
                    fixed_error_reduction_ratio = candidate_fixed_error / initial_fixed_error

                end


                ##########################################################################################
                #                            Compute Residuals                                           #
                ##########################################################################################
                
                # Begin "evaluate_unscaled_iteraiton_stats" #TODO: UPDATED! changed avg_primal_* and avg_dual_* to current_primal_* and current_dual_*
                
                while idx <= n_vars
                    original_primal_solution[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                    original_primal_gradient[LP, idx] = current_primal_gradient[LP, idx] * variable_rescaling[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                while idx <= current_LP_length
                    original_dual_solution[active_row + idx] = current_dual_solution[active_row + idx] / constraint_rescaling[active_row + idx]
                    original_primal_product[active_row + idx] = current_primal_product[active_row + idx] * constraint_rescaling[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x


                
                # Compute Convergence Information (will be used in Termination criteria check later on)
                # ++++++++++++++++++++++++++++++++++++++++++++++++ Primal Objective ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
                # Compute the primal objective
                while idx <= n_vars
                    shared_space[idx] = original_objective_vector[LP, idx] * original_primal_solution[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_primal_objective = original_objective_constant[LP] + shared_space[1]
                end
                # ++++++++++++++++++++++++++++++++++++++++++++++++ Primal Residual ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #                # Computing the primal residual (termination criteria 2)
                ## since variable bounds is part of the inequality constraints, we first calculate that residual
                ## Compute primal variable/constraint violations
                while idx <= n_vars
                    buffer_kkt_lower_variable_violation[LP, idx] = max(original_variable_lower_bounds[LP, idx] - original_primal_solution[LP, idx], 0.0)
                    buffer_kkt_upper_variable_violation[LP, idx] = max(original_primal_solution[LP, idx] - original_variable_upper_bounds[LP, idx], 0.0)
                    idx += block_stride
                end
                idx = threadIdx().x

                while idx <= n_vars
                    shared_space[idx] = abs(buffer_kkt_lower_variable_violation[LP, idx])^2 +
                                        abs(buffer_kkt_upper_variable_violation[LP, idx])^2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_l2_primal_residual = shared_space[1]
                end

                ## Calculating constraint violations (residuals)
                # based on this logic, a constraint is not violated if b - Ax < 0 because we have Ax > b
                while idx <= current_LP_length
                    shared_space[idx] = abs(max(original_right_hand_side[active_row + idx] - original_primal_product[active_row + idx], 0.0)) ^ 2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    CI_l2_primal_residual = sqrt(CI_l2_primal_residual + shared_space[1])
                end
                # ++++++++++++++++++++++++++++++++++++++++++++++ Dual Objective ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
                # Calculate the dual objective for termination criteria 1
                ## Compute dual variable as a dependant on primal and dual variables solution as r = c - K^T y
                while idx <= n_vars
                    buffer_kkt_reduced_costs[LP, idx] = max(original_primal_gradient[LP, idx], 0.0) * isfinite(original_variable_lower_bounds[LP, idx]) + 
                                                        min(original_primal_gradient[LP, idx], 0.0) * isfinite(original_variable_upper_bounds[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x

            
                ## Calculating the first term l^T r - u^T r
                while idx <= n_vars
                    if buffer_kkt_reduced_costs[LP, idx] > 0
                        shared_space[idx] = original_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                    elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                        shared_space[idx] = original_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                    else
                        shared_space[idx] = 0.0
                    end
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_dual_objective = original_objective_constant[LP] + shared_space[1]
                end

                ## Calculating the second term of dual objective function = q^T y
                while idx <= current_LP_length
                    shared_space[idx] = original_right_hand_side[active_row + idx] * original_dual_solution[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    CI_dual_objective += shared_space[1]
                end
                # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Dual Residual +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
                # Calculate the l2 dual residual (termination criteria 3)
                ## dual constraint violation check
                while idx <= n_vars
                    shared_space[idx] = abs(original_primal_gradient[LP, idx] - buffer_kkt_reduced_costs[LP, idx]) ^ 2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_l2_dual_residual = shared_space[1]
                end

                ## same logic as for primal, we calculate the dual variable bound violation (we need y>= 0) #TODO: ask why shadow price has to be >=0 in duality.
                while idx <= current_LP_length
                    shared_space[idx] = abs(max(-original_dual_solution[active_row + idx], 0.0))^2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    CI_l2_dual_residual = sqrt(CI_l2_dual_residual + shared_space[1])
                end

                
                # Compute Infeasibility Information (Ray-based infeasibility detection)
                
                sync_threads()
                
                while idx <= n_vars
                    shared_space[idx] = abs(original_primal_solution[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, var_stride, n_vars)
                # Everyone needs to know the primal ray norm
                primal_ray_norm = shared_space[1]

                # Calculate infeasibility primal solution and primal product
                if !iszero(primal_ray_norm)
                    while idx <= n_vars
                        buffer_kkt_primal_solution[LP, idx] = original_primal_solution[LP, idx] / primal_ray_norm
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    while idx <= current_LP_length
                        buffer_kkt_primal_product[active_row + idx] = original_primal_product[active_row + idx] / primal_ray_norm
                        idx += block_stride
                    end
                    idx = threadIdx().x
                else
                    while idx <= n_vars
                        buffer_kkt_primal_solution[LP, idx] = original_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    while idx <= current_LP_length
                        buffer_kkt_primal_product[active_row + idx] = original_primal_product[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                end

                # Compute infeasible variable/constraint violations
                while idx <= n_vars
                    if isfinite(original_variable_lower_bounds[LP, idx])
                        buffer_kkt_lower_variable_violation[LP, idx] = max(-buffer_kkt_primal_solution[LP, idx], 0.0)
                    else
                        buffer_kkt_lower_variable_violation[LP, idx] = 0.0
                    end
                    if isfinite(original_variable_upper_bounds[LP, idx])
                        buffer_kkt_upper_variable_violation[LP, idx] = max(buffer_kkt_primal_solution[LP, idx], 0.0)
                    else
                        buffer_kkt_upper_variable_violation[LP, idx] = 0.0
                    end
                    idx += block_stride
                end
                idx = threadIdx().x

                # Calculate the max primal ray infeasibility
                while idx <= n_vars
                    shared_space[idx] = max(abs(buffer_kkt_lower_variable_violation[LP, idx]), 
                                            abs(buffer_kkt_upper_variable_violation[LP, idx]))
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    II_max_primal_ray_infeasibility = shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = abs(max(-buffer_kkt_primal_product[active_row + idx], 0.0))
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    II_max_primal_ray_infeasibility = max(shared_space[1], II_max_primal_ray_infeasibility)
                end

                # Calculate the primal ray linear objective
                while idx <= n_vars
                    shared_space[idx] = original_objective_vector[LP, idx] * buffer_kkt_primal_solution[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    II_primal_ray_linear_objective = shared_space[1]
                end

                # Compute reduced costs and reduced costs violation
                while idx <= n_vars
                    buffer_kkt_reduced_costs[LP, idx] = max(original_primal_gradient[LP, idx] - original_objective_vector[LP, idx], 0.0) * isfinite(original_variable_lower_bounds[LP, idx]) + 
                                                        min(original_primal_gradient[LP, idx] - original_objective_vector[LP, idx], 0.0) * isfinite(original_variable_upper_bounds[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x

                # Compute the dual residual
                while idx <= n_vars
                    shared_space[idx] = abs(original_primal_gradient[LP, idx] - original_objective_vector[LP, idx] - buffer_kkt_reduced_costs[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    buffer_kkt_dual_res_inf = shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = abs(max(-original_dual_solution[active_row + idx], 0.0))
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    buffer_kkt_dual_res_inf = max(shared_space[1], buffer_kkt_dual_res_inf)
                end

                # Compute the dual objective
                while idx <= n_vars
                    if buffer_kkt_reduced_costs[LP, idx] > 0.0
                        shared_space[idx] = original_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                    elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                        shared_space[idx] = original_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                    else
                        shared_space[idx] = 0.0
                    end
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    buffer_kkt_dual_objective = shared_space[1] + original_objective_constant[LP]
                end
                while idx <= current_LP_length
                    shared_space[idx] = original_right_hand_side[active_row + idx] * original_dual_solution[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    buffer_kkt_dual_objective += shared_space[1]
                end

                # Compute infeasibility information using a scaling factor
                while idx <= n_vars
                    shared_space[idx] = abs(buffer_kkt_reduced_costs[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    scaling_factor = shared_space[1]
                end

                while idx <= current_LP_length
                    shared_space[idx] = abs(original_dual_solution[active_row + idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    scaling_factor = max(shared_space[1], scaling_factor)
                    if scaling_factor==0.0
                        II_max_dual_ray_infeasibility = 0.0
                        II_dual_ray_objective = 0.0
                    else
                        II_max_dual_ray_infeasibility = buffer_kkt_dual_res_inf / scaling_factor
                        II_dual_ray_objective = buffer_kkt_dual_objective / scaling_factor
                    end
                end

                total_iterations += Int32(termination_evaluation_frequency)

                ##########################################################################################
                #                            Check Termination Criteria                                  #
                ##########################################################################################
                sync_threads()

                # Check optimality criteria
                if idx==1
                    # CUDA.@cuprintln("sol = [$(original_primal_solution[1,1]), $(original_primal_solution[1,2])]")
                    # CUDA.@cuprintln("CI l2 dual residual = $(CI_l2_dual_residual * 1e6)")


                    # CUDA.@cuprintln("CI l2 primal residual = $(CI_l2_primal_residual * 1e6)")
                    # CUDA.@cuprintln("       - lower variable violation = $(buffer_kkt_lower_variable_violation[1,1]) , $(buffer_kkt_lower_variable_violation[1,2])")
                    # CUDA.@cuprintln("       - upper variable violation = $(buffer_kkt_upper_variable_violation[1,1]), $(buffer_kkt_upper_variable_violation[1,2])")

                    # CUDA.@cuprintln("       - lower_bound - sol = $(original_variable_lower_bounds[1,1] - original_primal_solution[1,1]), $(original_variable_lower_bounds[1,2] - original_primal_solution[1,2])")
                    # CUDA.@cuprintln("       - sol - upper_bound = $(original_primal_solution[1,1] - original_variable_upper_bounds[1,1]), $(original_primal_solution[1,2] - original_variable_upper_bounds[1,2])")

                    # # original_right_hand_side[active_row + idx] - original_primal_product[active_row + idx]

                    # CUDA.@cuprintln("       - b - Ax [1] = $(original_right_hand_side[1] - original_primal_product[1])")
                    # CUDA.@cuprintln("       - b - Ax [2] = $(original_right_hand_side[2] - original_primal_product[2])")
                    # CUDA.@cuprintln("       - b - Ax [3] = $(original_right_hand_side[3] - original_primal_product[3])")
                    # CUDA.@cuprintln("       - b - Ax [4] = $(original_right_hand_side[4] - original_primal_product[4])")
                    # CUDA.@cuprintln("       - b - Ax [5] = $(original_right_hand_side[5] - original_primal_product[5])")
                    # CUDA.@cuprintln("       - b - Ax [6] = $(original_right_hand_side[6] - original_primal_product[6])")

                    

                    # CUDA.@cuprintln("CI l2 primal objective = $(CI_primal_objective)")
                    # CUDA.@cuprintln("       - c[1] * sol[1] = $(original_objective_vector[1,1] * original_primal_solution[1,1])")
                    # CUDA.@cuprintln("       - c[2] * sol[2] = $(original_objective_vector[1,2] * original_primal_solution[1,2])")


                    # CUDA.@cuprintln("CI l2 dual objective = $(CI_dual_objective)")
                    # Check if the dual objective value is above the B&B global upper bound and we
                    # satisfy the tolerance for dual feasibility
                    # (if we're past the first 10 iterations)
                    if (CI_dual_objective > global_upper_bound + abs_tol) &&
                        (CI_l2_dual_residual < abs_tol + rel_tol*cache_l2_norm_primal_linear_objective) 
                        termination_reason[LP] = TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT
                    end

                    # If we want to skip harder-than-average problems, and we've already solved at least
                    # 100 problems, check if the current number of iterations is over two times the running
                    # average
                    if skip_hard_problems && (unsafe_load(CUDA.pointer(global_counter, 1)) > 100)
                        if total_iterations > 2*(unsafe_load(CUDA.pointer(iteration_counter, 1))/unsafe_load(CUDA.pointer(global_counter, 1)))
                            termination_reason[LP] = TERMINATION_REASON_IMPATIENCE
                        end
                    end

                    # Check iteration limit
                    if total_iterations >= iteration_limit
                        termination_reason[LP] = TERMINATION_REASON_ITERATION_LIMIT
                    end

                    # Check KKT matrix pass limit
                    # if cumulative_kkt_passes >= kkt_matrix_pass_limit
                    #     termination_reason[LP] = TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT
                    # end

                    # Check for numerical errors
                    if numerical_error[1]
                        termination_reason[LP] = TERMINATION_REASON_NUMERICAL_ERROR
                    end

                    # Check if we're within the tolerances for primal and dual infeasibility, and that there's
                    # a sufficiently small gap between the primal and dual objective values.

                    if (CI_l2_dual_residual < abs_tol + rel_tol*cache_l2_norm_primal_linear_objective) &&
                        (CI_l2_primal_residual < abs_tol + rel_tol*cache_l2_norm_primal_right_hand_side) &&
                        (abs(CI_primal_objective - CI_dual_objective) < abs_tol + rel_tol*(abs(CI_primal_objective)+abs(CI_dual_objective))) 

                        termination_reason[LP] = TERMINATION_REASON_OPTIMAL
                    end
                    
                    # Checking on infeasibility after convergence failure
                    ## Check primal infeasibility (if we're past the first 10 iterations)
                    if (II_dual_ray_objective > 0.0) &&
                        ((II_max_dual_ray_infeasibility / II_dual_ray_objective) <= eps_primal_infeasible) && (total_iterations > 10)

                        termination_reason[LP] = TERMINATION_REASON_PRIMAL_INFEASIBLE
                    end

                    # Check dual infeasibility (if we're past the first 10 iterations)
                    if (II_primal_ray_linear_objective < 0.0) && 
                        ((II_max_primal_ray_infeasibility / (-II_primal_ray_linear_objective)) <= eps_dual_infeasible)  && (total_iterations > 10)

                        termination_reason[LP] = TERMINATION_REASON_DUAL_INFEASIBLE
                    end
                end

                
                # Update Solutions (UPDATED!)
                

                # If the LP is finished, update the solutions and objective(s), and
                # then break out of the iteration loop and move on to the next LP
                sync_threads()

                reason = unsafe_load(CUDA.pointer(termination_reason, LP))
                if reason == TERMINATION_REASON_OPTIMAL
                    while idx <= n_vars 
                        solutions[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        if return_code==Int32(1)
                            objectives[LP] = CI_primal_objective
                        elseif return_code==Int32(2)
                            objectives[LP] = CI_dual_objective
                        else
                            objectives[LP, Int32(1)] = CI_primal_objective
                            objectives[LP, Int32(2)] = CI_dual_objective
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_PRIMAL_INFEASIBLE) || 
                    (reason == TERMINATION_REASON_DUAL_INFEASIBLE)
                    while idx <= n_vars 
                        solutions[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        if return_code != Int32(3)
                            objectives[LP] = Inf
                        else
                            objectives[LP, Int32(1)] = Inf
                            objectives[LP, Int32(2)] = Inf
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_ITERATION_LIMIT) || 
                    (reason == TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT) ||
                    (reason == TERMINATION_REASON_NUMERICAL_ERROR)
                    while idx <= n_vars 
                        solutions[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        if return_code != Int32(3)
                            objectives[LP] = -Inf
                        else
                            objectives[LP, Int32(1)] = -Inf
                            objectives[LP, Int32(2)] = -Inf
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT)
                    while idx <= n_vars
                        solutions[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        # Special case: This checks for a feasible dual objective value
                        # above the global upper bound in a B&B algorithm. The primal
                        # objective value may not be meaningful/feasible.
                        if return_code==Int32(1)
                            objectives[LP] = CI_primal_objective
                        elseif return_code==Int32(2)
                            objectives[LP] = CI_dual_objective
                        else
                            objectives[LP, Int32(1)] = CI_primal_objective
                            objectives[LP, Int32(2)] = CI_dual_objective
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_IMPATIENCE)
                    while idx < n_vars
                        solutions[LP, idx] = current_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        if return_code != Int32(3)
                            objectives[LP] = -Inf
                        else
                            objectives[LP, Int32(1)] = -Inf
                            objectives[LP, Int32(2)] = -Inf
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    LP += grid_stride
                    break
                end

                ##########################################################################################
                #                            Check Restart Conditions                                    #
                ##########################################################################################

                # Check if we need to do a restart
                if idx == 1
                    
                    if fixed_error_reduction_ratio < necessary_reduction_for_restart

                        if fixed_error_reduction_ratio < sufficient_reduction_for_restart
                            # CUDA.@cuprintln("We are doing restart because of sufficient reduction")
                            do_restart[1] = true
                        elseif fixed_error_reduction_ratio > last_reduction_ratio
                            # CUDA.@cuprintln("We are doing restart because of last reduction ratio")
                            do_restart[1] = true
                        end

                        last_reduction_ratio = fixed_error_reduction_ratio

                    elseif inner_iterations >= artificial_ratio_for_restart * (total_iterations - 1)
                        # CUDA.@cuprintln("We are doing restart because of long inner loop")
                        do_restart[1] = true
                    end
                end
                ##########################################################################################
                #                            Perform Restart if do_restart = true                        #
                ##########################################################################################
                sync_threads()
                
                if unsafe_load(CUDA.pointer(do_restart, 1)) == true 
                    
                    # if idx == 1
                    #     termination_reason[LP] = TERMINATION_REASON_NUMERICAL_ERROR
                    # end
                    # break
                    # Primal distance
                    while idx <= n_vars 
                        shared_space[idx] = (current_primal_solution[LP, idx] - initial_primal_solution[LP, idx]) ^ 2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)

                    if idx == 1
                        restart_primal_distance = sqrt(primal_weight[1]) * sqrt(shared_space[1])
                    end
                    
                    # Dual distance
                    while idx <= current_LP_length
                        shared_space[idx] = (current_dual_solution[active_row + idx] - initial_dual_solution[active_row + idx]) ^ 2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)

                    if idx == 1
                        restart_dual_distance = (1 / sqrt(primal_weight[1])) * sqrt(shared_space[1])
                        restart_error = log(restart_primal_distance / restart_dual_distance )
                        CUDA.@cuprintln("dual dist = $(restart_dual_distance), primal dist = $(restart_primal_distance), primal weight = $(primal_weight[1]), error = $restart_error")
                        # restart_error = log(restart_dual_distance / restart_primal_distance)
                        # old_error = log(restart_dual_distance / restart_primal_distance)
                        sum_restart_error += restart_error
                    end


                    #### Update Information About the Last Restart
                    ## resetting step itertation back to 0 
                    inner_iterations = Int32(0)

                    # initialize the primal anchor 
                    while idx <= n_vars
                        initial_primal_solution[LP, idx] = current_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # initialize the dual anchor
                    while idx <= current_LP_length
                        initial_dual_solution[active_row + idx] = current_dual_solution[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # # # updating the current_primal_solution and current_dual_solution 
                    # while idx <= n_vars 
                    #     current_primal_solution[LP, idx] = pdhg_primal_solution[LP, idx]
                    #     idx += block_stride
                    # end
                    # idx = threadIdx().x

                    # # initialize the dual anchor
                    # while idx <= current_LP_length
                    #     current_dual_solution[active_row + idx] = pdhg_dual_solution[active_row + idx]
                    #     idx += block_stride
                    # end
                    # idx = threadIdx().x

                    # we do not restart to average in Halpern Reflection Scheme
                    if idx==1
                        restart_choice = RESTART_CHOICE_LAST_ITERATE_RESET
                    end
                else
                    if idx==1
                        restart_choice = RESTART_CHOICE_NO_RESTART
                    end
                end
            


                ###################################################################
                ##### Compute a New Primal Weight (if a restart was used) 
                ###################################################################
                if idx==1
                    if restart_choice == RESTART_CHOICE_LAST_ITERATE_RESET
                        if restart_primal_distance > eps() && restart_dual_distance > eps()
                            # CUDA.@cuprintln("initial weight = $(primal_weight[1])")
                            
                            primal_weight[1] = exp(log(primal_weight[1]) - pid_KP * restart_error - pid_KI * sum_restart_error - pid_KD * (restart_error - last_restart_error))
                            
                            # primal_weight[1] = 0.1
                            # primal_weight[1] = exp(0.5 * old_error + (1 - 0.5) * log(primal_weight[1]))
                            # primal_weight[1] *= exp(pid_KP * restart_error + pid_KI * sum_restart_error + pid_KD * (restart_error - last_restart_error))
                            # CUDA.@cuprintln("Restart : Updated primal weight = $(primal_weight[1]) using the error $restart_error, and sum of error = $(sum_restart_error) and error diff = $(restart_error - last_restart_error)")
                        end
                    end
                    # update last restart error 
                    last_restart_error = restart_error
                end
                sync_threads()
        end
    end

    return nothing
end

    #-----------------------------------------------------------------------------------------

    # # Add to the cumulative kkt pass count
    # if idx==1
    #     cumulative_kkt_passes += 2.0
    # end

# A quick function to calculate a parallel sum over the first max_len elements in the shared space
function parallel_sum(shared, block_stride, reduction_stride, maxlen)
    sync_threads()
    idx = threadIdx().x
    while reduction_stride > 0
        while idx <= reduction_stride && idx + reduction_stride <= maxlen
            shared[idx] += shared[idx + reduction_stride]
            idx += block_stride
        end
        idx = threadIdx().x
        sync_threads()
        reduction_stride >>>= 1
    end
    return nothing
end

# A quick function to calculate a parallel max over the first max_len elements in the shared space
function parallel_max(shared, block_stride, reduction_stride, maxlen)
    sync_threads()
    idx = threadIdx().x
    while reduction_stride > 0
        while idx <= reduction_stride && idx + reduction_stride <= maxlen
            shared[idx] = max(shared[idx], shared[idx + reduction_stride])
            idx += block_stride
        end
        idx = threadIdx().x
        sync_threads()
        reduction_stride >>>= 1
    end
    return nothing
end
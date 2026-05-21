
function main_loop_kernel(
    solutions,                        # [n_LPs, n_vars−1] = Holds solutions to LPs #TODO: I need to check this, I think it should be [n_LPs, n_vars]
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
    scaled_objective_constant,        # [n_LPs] = Objective constant after initial rescaling
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
    buffer_primal_gradient,           # [n_LPs, n_vars] = 
    initial_primal_solution,
    initial_dual_solution,
    next_primal_solution,
    next_dual_solution,
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
    delta_primal_halpern,
    delta_dual_halpern,
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
        end

        # Other information that every thread needs (stored in the Registery memory)
            ## The values are always identical across threads — it is replicated rather than shared for performance purposes.
        step_iterations = Int32(1) # inner loop counter
        iteration = Int32(1) # total iteration counter (T)

        
        primal_ray_norm = 0.0


        
        # Information that is much easier to save as static shared memory (L1 Block memory)
        do_restart = @cuStaticSharedMem(Bool, 1)
        primal_weight = @cuStaticSharedMem(Float64, 1)
        numerical_error = @cuStaticSharedMem(Bool, 1)
        if idx==1
            do_restart[1] = false
            primal_weight[1] = input_primal_weight[LP]
            step_size[1] = input_step_size[LP]
            numerical_error[1] = false
        end
        
        # Set up the starting row for this LP (minus 1, so that the first
        # row to consider is `active_row + 1`)
        active_row = (LP-Int32(1)) * total_LP_length 

        # Set up current values to match inputs
        while idx <= n_vars
            buffer_primal_gradient[LP, idx] = scaled_objective_vector[LP, idx]
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
        while true

            ###################################################################
            ##### Phase 1 of the Main Loop : Check termination at each iteration 
            ###################################################################
    
            
            # Add to the cumulative kkt pass count
            if idx==1
                cumulative_kkt_passes += 2.0
            end

            # Update the average primal/dual solutions, average primal gradient,
            # and average primal product
            #TODO: UPDATED! Completely removed since no averaging is happening in cuPDLPx

            sync_threads()
            

            ###################################################################
            ##### Begin "evaluate_unscaled_iteraiton_stats" #TODO: UPDATED! changed avg_primal_* and avg_dual_* to current_primal_* and current_dual_*
            ###################################################################
            while idx <= n_vars
                original_primal_solution[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                original_primal_gradient[LP, idx] = buffer_primal_gradient[LP, idx] * variable_rescaling[LP, idx]
                idx += block_stride
            end
            idx = threadIdx().x
            while idx <= current_LP_length
                original_dual_solution[active_row + idx] = current_dual_solution[active_row + idx] / constraint_rescaling[active_row + idx]
                original_primal_product[active_row + idx] = current_primal_product[active_row + idx] * constraint_rescaling[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x


            ###################################################################
            ##### Compute Convergence Information
            ###################################################################

            # Compute primal variable/constraint violations
            while idx <= n_vars
                buffer_kkt_lower_variable_violation[LP, idx] = max(original_variable_lower_bounds[LP, idx] - original_primal_solution[LP, idx], 0.0)
                buffer_kkt_upper_variable_violation[LP, idx] = max(original_primal_solution[LP, idx] - original_variable_upper_bounds[LP, idx], 0.0)
                idx += block_stride
            end
            idx = threadIdx().x

            # Compute the primal objective and residual
            while idx <= n_vars
                shared_space[idx] = original_objective_vector[LP, idx] * original_primal_solution[LP, idx]
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                CI_primal_objective = original_objective_constant[LP] + shared_space[1]
            end

            # CI_l_inf_primal_residual
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
            while idx <= current_LP_length
                shared_space[idx] = abs(max(original_right_hand_side[active_row + idx] - original_primal_product[active_row + idx], 0.0))^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                CI_l2_primal_residual = sqrt(CI_l2_primal_residual + shared_space[1])
            end
            
            # Compute dual stats
            while idx <= n_vars
                buffer_kkt_reduced_costs[LP, idx] = max(original_primal_gradient[LP, idx], 0.0) * isfinite(original_variable_lower_bounds[LP, idx]) + 
                                                    min(original_primal_gradient[LP, idx], 0.0) * isfinite(original_variable_upper_bounds[LP, idx])
                idx += block_stride
            end
            idx = threadIdx().x

            # Calculate the dual objective
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
            while idx <= current_LP_length
                shared_space[idx] = original_right_hand_side[active_row + idx] * original_dual_solution[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                CI_dual_objective += shared_space[1]
            end

            # Calculate the l2 dual residual
            while idx <= n_vars
                shared_space[idx] = abs(original_primal_gradient[LP, idx] - buffer_kkt_reduced_costs[LP, idx])^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                CI_l2_dual_residual = shared_space[1]
            end
            while idx <= current_LP_length
                shared_space[idx] = abs(max(-original_dual_solution[active_row + idx], 0.0))^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                CI_l2_dual_residual = sqrt(CI_l2_dual_residual + shared_space[1])
            end



            ###################################################################
            ##### Compute Infeasibility Information 
            ###################################################################

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


            ###################################################################
            ##### Check Termination Criteria
            ###################################################################

            # Compute the current primal residual
            while idx <= n_vars
                buffer_kkt_lower_variable_violation[LP, idx] = max(scaled_variable_lower_bounds[LP, idx] - current_primal_solution[LP, idx], 0.0)
                buffer_kkt_upper_variable_violation[LP, idx] = max(current_primal_solution[LP, idx] - scaled_variable_upper_bounds[LP, idx], 0.0)
                idx += block_stride
            end
            idx = threadIdx().x

            ## Compute the current primal objective and l2 primal residual
            while idx <= n_vars
                shared_space[idx] = scaled_objective_vector[LP, idx] * current_primal_solution[LP, idx]
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                primal_objective_storage = shared_space[1] + scaled_objective_constant[LP]
            end

            while idx <= n_vars
                shared_space[idx] = abs(buffer_kkt_lower_variable_violation[LP, idx])^2 +
                                    abs(buffer_kkt_upper_variable_violation[LP, idx])^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                l2_primal_residual = shared_space[1]
            end
            while idx <= current_LP_length
                shared_space[idx] = abs(max(scaled_right_hand_side[active_row + idx] - current_primal_product[active_row + idx], 0.0))^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                l2_primal_residual = sqrt(shared_space[1] + l2_primal_residual)
            end

            ## Compute the current dual residual
            while idx <= n_vars
                buffer_kkt_reduced_costs[LP, idx] = max(buffer_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_lower_bounds[LP, idx]) + 
                                                    min(buffer_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_upper_bounds[LP, idx])
                idx += block_stride
            end
            idx = threadIdx().x

            ## Calculate the current dual objective
            while idx <= n_vars
                if buffer_kkt_reduced_costs[LP, idx] > 0.0
                    shared_space[idx] = scaled_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                    shared_space[idx] = scaled_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                else
                    shared_space[idx] = 0.0
                end
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                buffer_kkt_dual_objective = scaled_objective_constant[LP] + shared_space[1]
            end
            while idx <= current_LP_length
                shared_space[idx] = scaled_right_hand_side[active_row + idx] * current_dual_solution[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                buffer_kkt_dual_objective += shared_space[1]
            end
            
            # Calculate the l2 dual residual
            while idx <= n_vars
                shared_space[idx] = abs(buffer_primal_gradient[LP, idx] - buffer_kkt_reduced_costs[LP, idx])^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                l2_dual_residual = shared_space[1]
            end
            while idx <= current_LP_length
                shared_space[idx] = abs(max(-current_dual_solution[active_row + idx], 0.0))^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                l2_dual_residual += shared_space[1]
            end

            # Calculate the "current" kkt residual
            if idx==1
                current_kkt_residual = sqrt(primal_weight[1] * l2_primal_residual^2 + 
                                            (1/primal_weight[1]) * l2_dual_residual^2 + 
                                            abs(primal_objective_storage - buffer_kkt_dual_objective)^2)
            end

            # Compute the current primal residual
            while idx <= n_vars
                buffer_kkt_lower_variable_violation[LP, idx] = max(scaled_variable_lower_bounds[LP, idx] - avg_primal_solution[LP, idx], 0.0)
                buffer_kkt_upper_variable_violation[LP, idx] = max(avg_primal_solution[LP, idx] - scaled_variable_upper_bounds[LP, idx], 0.0)
                idx += block_stride
            end
            idx = threadIdx().x

            # Compute the current primal objective and l2 primal residual
            while idx <= n_vars
                shared_space[idx] = scaled_objective_vector[LP, idx] * avg_primal_solution[LP, idx]
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                primal_objective_storage = shared_space[1] + scaled_objective_constant[LP]
            end

            while idx <= n_vars
                shared_space[idx] = abs(buffer_kkt_lower_variable_violation[LP, idx])^2 + 
                                    abs(buffer_kkt_upper_variable_violation[LP, idx])^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                l2_primal_residual = shared_space[1]
            end
            while idx <= current_LP_length
                shared_space[idx] = abs(max(scaled_right_hand_side[active_row + idx] - avg_primal_product[active_row + idx], 0.0))^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                l2_primal_residual = sqrt(shared_space[1] + l2_primal_residual)
            end

            # Compute the current dual residual
            while idx <= n_vars
                buffer_kkt_reduced_costs[LP, idx] = max(avg_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_lower_bounds[LP, idx]) + 
                                                    min(avg_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_upper_bounds[LP, idx])
                idx += block_stride
            end
            idx = threadIdx().x

            # Calculate the current dual objective
            while idx <= n_vars
                if buffer_kkt_reduced_costs[LP, idx] > 0.0
                    shared_space[idx] = scaled_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                    shared_space[idx] = scaled_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                else
                    shared_space[idx] = 0.0
                end
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                buffer_kkt_dual_objective = scaled_objective_constant[LP] + shared_space[1]
            end
            while idx <= current_LP_length
                shared_space[idx] = scaled_right_hand_side[active_row + idx] * avg_dual_solution[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                buffer_kkt_dual_objective += shared_space[1]
            end

            # Calculate the l2 dual residual
            while idx <= n_vars
                shared_space[idx] = abs(avg_primal_gradient[LP, idx] - buffer_kkt_reduced_costs[LP, idx])^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx==1
                l2_dual_residual = shared_space[1]
            end
            while idx <= current_LP_length
                shared_space[idx] = abs(max(-avg_dual_solution[active_row + idx], 0.0))^2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                l2_dual_residual += shared_space[1]
            end
            sync_threads()

            # Check optimality criteria
            if idx==1
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
                    if iteration > 2*(unsafe_load(CUDA.pointer(iteration_counter, 1))/unsafe_load(CUDA.pointer(global_counter, 1)))
                        termination_reason[LP] = TERMINATION_REASON_IMPATIENCE
                    end
                end

                # Check iteration limit
                if iteration >= iteration_limit
                    termination_reason[LP] = TERMINATION_REASON_ITERATION_LIMIT
                end

                # Check KKT matrix pass limit
                if cumulative_kkt_passes >= kkt_matrix_pass_limit
                    termination_reason[LP] = TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT
                end

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
                
                # Check primal infeasibility (if we're past the first 10 iterations)
                if (II_dual_ray_objective > 0.0) &&
                    ((II_max_dual_ray_infeasibility / II_dual_ray_objective) <= eps_primal_infeasible) &&
                    (iteration >= 10)
                    termination_reason[LP] = TERMINATION_REASON_PRIMAL_INFEASIBLE
                end

                # Check dual infeasibility (if we're past the first 10 iterations)
                if (II_primal_ray_linear_objective < 0.0) && 
                    ((II_max_primal_ray_infeasibility / (-II_primal_ray_linear_objective)) <= eps_dual_infeasible) &&
                    (iteration >= 10)
                    termination_reason[LP] = TERMINATION_REASON_DUAL_INFEASIBLE
                end
            end


            ###################################################################
            ##### Update Solutions (UPDATED!)
            ###################################################################

            # If the LP is finished, update the solutions and objective(s), and
            # then break out of the iteration loop and move on to the next LP
            sync_threads()
            reason = unsafe_load(CUDA.pointer(termination_reason, LP))
            if reason == TERMINATION_REASON_OPTIMAL
                while idx < n_vars # Intentionally skipping the epigraph variable
                    solutions[LP, idx] = current_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
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
                    iterations[LP] = iteration
                    CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                    CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                end
                LP += grid_stride
                break
            elseif (reason == TERMINATION_REASON_PRIMAL_INFEASIBLE) || 
                (reason == TERMINATION_REASON_DUAL_INFEASIBLE)
                while idx < n_vars # Intentionally skipping the epigraph variable
                    solutions[LP, idx] = current_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
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
                    iterations[LP] = iteration
                    CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                    CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                end
                LP += grid_stride
                break
            elseif (reason == TERMINATION_REASON_ITERATION_LIMIT) || 
                (reason == TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT) ||
                (reason == TERMINATION_REASON_NUMERICAL_ERROR)
                while idx < n_vars # Intentionally skipping the epigraph variable
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
                    iterations[LP] = iteration
                    CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                    CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                end
                LP += grid_stride
                break
            elseif (reason == TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT)
                while idx < n_vars
                    solutions[LP, idx] = current_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
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
                    iterations[LP] = iteration
                    CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                    CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
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
                    iterations[LP] = iteration
                    CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                    CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                end
                LP += grid_stride
                break
            end
            

            ###################################################################
            ##### Update Buffer Primal Gradient 
            ###################################################################
            while idx <= n_vars
                buffer_primal_gradient[LP, idx] = scaled_objective_vector[LP, idx] - current_dual_product[LP, idx]
                idx += block_stride
            end
            idx = threadIdx().x


            ###################################################################
            ##### Check for restart criteria UPDATED! fixed-point error method has been implemented!
            ###################################################################
            sync_threads()

            if idx==1
                # Initialize do_restart
                do_restart[1] = false

                if step_iterations >= artificial_ratio_for_restart * (iteration - 1)
                    do_restart[1] = true
                end
            end

                
            if unsafe_load(CUDA.pointer(do_restart, 1))==false
            
                if idx==1

                    initial_fixed_error = sqrt((primal_weight / step_size) * anchor_squared_delta_primal + (1 / (step_size * primal_weight)) * anchor_squared_delta_dual + 2 * anchor_cross_term)
                    candidate_fixed_error = sqrt((primal_weight / step_size) * squared_delta_primal + (1 / (step_size * primal_weight)) * squared_delta_dual + 2 * cross_term)
                    fixed_error_reduction_ratio = candidate_fixed_error / initial_fixed_error

                    # Check if we need to do a restart
                    if fixed_error_reduction_ratio < necessary_reduction_for_restart
                        if fixed_error_reduction_ratio < sufficient_reduction_for_restart
                            do_restart[1] = true
                        elseif fixed_error_reduction_ratio > last_reduction_ratio
                            do_restart[1] = true
                        end
                    end
                    
                    last_reduction_ratio = fixed_error_reduction_ratio
                end
                
            end


            #>>>###################################################################
            #>>>##### Apply the Reset (if do_restart==true) #TODO: UPDATED! step_iterations will reset to one, the anchor (initial_* will get updated), current solutions will be set to T(current_*) (smae values as the anchor)
            #>>>###################################################################
            sync_threads()

            if unsafe_load(CUDA.pointer(do_restart, 1)) == true 

                # Primal distance
                while idx <= n_vars 
                    shared_space[idx] = (current_primal_solution[LP, idx] - initial_primal_solution[LP, idx]) ^ 2
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)

                if idx == 1
                        restart_primal_distance = sqrt(primal_weight[1]) * sqrt(shared_space[1])
                end
                
                # Dual distance

                while idx <= current_LP_length
                    shared_space[idx] = (current_dual_solution[LP, idx] - initial_dual_solution[LP, idx]) ^ 2
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)

                if idx == 1
                    restart_dual_distance = (1 / sqrt(primal_weight[1])) * sqrt(shared_space[1])
                end

                # Calculating the error based on distance 
                if idx == 1
                    restart_error = log(restart_primal_distance / restart_dual_distance )
                end


                # udpating the sum 
                if idx == 1
                    sum_restart_error += restart_error
                end


                #### Update Information About the Last Restart
                ## resetting step itertation back to one 
                step_iterations = Int32(1)

                # initialize the primal anchor 
                while idx <= n_vars
                    initial_primal_solution[LP, idx] = next_primal_solution[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x

                # initialize the dual anchor
                while idx <= current_LP_length
                    initial_dual_solution[active_row + idx] = next_dual_solution[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x

                # updating the current_primal_solution and current_dual_solution 
                while idx <= n_vars 
                    current_primal_solution[LP, idx] = next_primal_solution[LP, idx]
                end
                idx = threadIdx().x

                # initialize the dual anchor
                while idx <= current_LP_length
                    current_dual_solution[active_row + idx] = next_dual_solution[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x

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
            ##### Compute a New Primal Weight (if a restart was used) #TODO: UPDATED! to have pid controller 
            ###################################################################
            if idx==1
                if restart_choice != RESTART_CHOICE_NO_RESTART
                    
                    primal_weight[1] = log(primal_weight[1]) - pid_KP * restart_error - pid_KI * sum_restart_error - pid_KD * (restart_error - last_restart_error)
                end
                # update last restart error 
                last_restart_error = restart_error
            end



            #########################################################################
            ##### Phase 2 of the Main Loop : Taking one Halpern Reflected PDHG Step #
            #########################################################################

            sync_threads()
            local_primal_weight = unsafe_load(CUDA.pointer(primal_weight, 1))
            local_step_size = unsafe_load(CUDA.pointer(step_size, 1))
            
            # Step 1) Add one to the total iterations tracker and cumulative kkt passes
            
            if idx==1
                cumulative_kkt_passes += Int32(1)
            end
                
            # Step 2) Operate on matrices to:
            # A) Update next_primal_solution (next_x= T(x[k])) -> next_primal_solution will be used to initialze the anchor if restart is triggered

            while idx <= n_vars 
                next_primal_solution[LP, idx] = min(scaled_variable_upper_bounds[LP, idx], max(scaled_variable_lower_bounds[LP, idx], current_primal_solution[LP, idx] - 
                                                (local_step_size/local_primal_weight) * (scaled_objective_vector[LP, idx] - current_dual_product[LP, idx])))
            end

            idx = threadIdx().x

            # Compute delta primal for dual step
            while idx <= n_vars 
                delta_primal[LP, idx] = next_primal_solution[LP, idx] - current_primal_solution[LP, idx]
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

            # Compute other terms normally
            while idx <= current_LP_length
                delta_dual[active_row + idx] = max(0.0, (current_dual_solution[active_row + idx] + 
                                    (local_primal_weight * local_step_size) * (scaled_right_hand_side[active_row + idx] -
                                    ((Int32(1) + extrapolation_coefficient) * 
                                    delta_primal_product[active_row + idx]) - 
                                    (extrapolation_coefficient * current_primal_product[active_row + idx])))) - 
                                    current_dual_solution[active_row + idx]

                shared_space[idx] = delta_primal_product[active_row + idx] * delta_dual[active_row + idx]
                idx += block_stride
            end
            
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            if idx==1
                cross_term = shared_space[1]
                if step_iterations == 1
                    anchor_cross_term  = cross_term
                end
            end
            # Getting next dual solution using delta_duals next_y = T(y) -> next_dual_solution will be used to initialze the anchor if restart is triggered
            while idx <= current_LP_length 
                next_dual_solution[active_row + idx] = delta_dual[active_row + idx] + current_dual_solution[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x
            
            # taking the Halpern scheme with Reflection step
            while idx <= n_vars
                current_primal_solution[LP, idx] = ((step_iterations + 1)/ (step_iterations + 2)) * ( ( 1 + reflection_coefficient) * next_primal_solution[LP, idx] - reflection_coefficient * current_primal_solution[LP, idx]) + (1 / (step_iterations + 2)) * initial_primal_solution[LP, idx]
                idx += block_stride
            end
            idx = threadIdx().x

            

            # Step 4) Update the solution
            sync_threads()

            
            while idx <= current_LP_length
                current_primal_product[active_row + idx] += delta_primal_product[active_row + idx]
                current_dual_solution[active_row + idx] = ((step_iterations + 1)/ (step_iterations + 2)) * ( ( 1 + reflection_coefficient) * next_dual_solution[active_row + idx] - reflection_coefficient * current_dual_solution[active_row + idx]) + (1 / (step_iterations + 2)) * initial_dual_solution[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x

            
            # Compute delta primal and delta dual Δz = (Δx, Δy) with Δz = HalpernReflected(z) - PDHG(z) for fixed-point error calculation
            while idx <= n_vars 
                delta_primal_halpern[LP, idx] = current_primal_solution[LP, idx] - next_primal_solution[LP, idx]
                idx += block_stride
            end

            idx = threadIdx().x


            # Compute the l2 norm of delta_primal and delta_dual for CURRENT solution and if we are the at first inner iteration, storing them for anchor
            while idx <= n_vars 
                shared_space[idx] = delta_primal_halpern[LP, idx] ^ 2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, var_stride, n_vars)
            if idx == 1
                squared_delta_primal = shared_space[1]
                if step_iterations == 1
                    anchor_squared_delta_primal = squared_delta_primal
                end
            end

            while idx <= current_LP_length
                delta_dual_halpern[active_row + idx] = current_dual_solution[active_row + idx] - next_dual_solution[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x

            while idx <= current_LP_length
                shared_space[idx] = delta_dual_halpern[active_row + idx] ^ 2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)

            if idx == 1
                squared_delta_dual = shared_space[1]
                if step_iterations == 1
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

            step_iterations += 1
            iteration += 1
        end
    end
    return nothing
end

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
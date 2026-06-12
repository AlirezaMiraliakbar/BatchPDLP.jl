
function PDLP(
    PDLP_data::PDLPData; 
    solutions::CuArray{Float64}=CuArray{Float64}(undef, PDLP_data.dims.n_LPs, PDLP_data.dims.n_vars), 
    objectives::CuArray{Float64}=CuArray{Float64}(undef, PDLP_data.dims.n_LPs),
    return_dual_obj::Bool=false,
    return_both_obj::Bool=false,
    global_upper_bound::Float64=Inf
    )

    # Check return conditions and verify that `objectives` storage is correctly sized
    if return_dual_obj + return_both_obj == 2
        error("Only one return condition allowed")
    end

    if !return_dual_obj && !return_both_obj
        return_code = Int32(1)
        if size(objectives, 2) != 1
            error("Objective storage sized incorrectly")
        end
    elseif return_dual_obj
        return_code = Int32(2)
        if size(objectives, 2) != 1
            error("Objective storage sized incorrectly")
        end
    elseif return_both_obj
        return_code = Int32(3)
        if size(objectives, 2) != 2
            error("Objective storage sized incorrectly")
        end
    end
    
    # # println("total LP length = $(PDLP_data.dims.total_LP_length)")
    # LP = 576
    # println("constriant matrix LP $LP: \n")
    # total_LP_length = PDLP_data.dims.total_LP_length
    # LP_section = 1 + (LP - 1) * total_LP_length : LP * total_LP_length
    # println(PDLP_data.original_problem.constraint_matrix[LP_section,:])

    # println("variable lower bounds of $LP: \n")
    # println(PDLP_data.original_problem.variable_lower_bounds[LP, :])

    # println("variable upper bounds of $LP: \n")
    # println(PDLP_data.original_problem.variable_upper_bounds[LP, :])
    # # println(size(PDLP_data.original_problem.variable_upper_bounds))
    # println("right hand side $LP: \n")
    # println(PDLP_data.original_problem.right_hand_side[LP_section])

    # println("objective vector $LP: \n")
    # println(PDLP_data.original_problem.objective_vector[LP, :])

    # println("objective constant $LP: \n")
    # println(PDLP_data.original_problem.objective_constant[LP, :])
    # error()   

    # Validate the LP data we've been given to make sure the numbers are all valid and the dimensions
    # of participating matrices are correct
    println("Validating the PDLP_data...")
    validate(PDLP_data)

    # Reset all fields relevant to problem status
    println("Resetting all fields...")
    reset_all_fields!(PDLP_data)

    # Perform rescaling. Note that if hot-starting is to be added in the future, you should save all
    # primal/dual information (rather than only primal solutions and dual objectives), and then here,
    # instead of rescaling the problem immediately, first un-scale the primal/dual statuses of the
    # previous run, then do re-scaling, then scale the primal/dual statuses again. This will also
    # require adding storage for these values, and un-commenting some lines in `main_loop.jl`
    # to allow hot-starting to impact the main PDLP algorithm. 

    # rescaling happens of original_problem -> ruiz_scaling (ruiz_var_kernel, ruiz_const_kernel, scale_problem) -> pock_chambolle -> scaled_problem
    # TODO: our constraint scaling stuff at ruiz_scaling happens inside scaling_part2_kernel that we are interesed in...
    println("Scaling the problem...")
    rescale_problem(
        PDLP_data.original_problem, 
        PDLP_data.scaled_problem, 
        PDLP_data.variable_rescaling, 
        PDLP_data.constraint_rescaling, 
        PDLP_data.dims,
        PDLP_data.parameters
        )
    # data from cuPDLPx for subproblem 18 for debugging
    # CUDA.@allowscalar begin
    #     PDLP_data.scaled_problem.constraint_matrix[1,1] = 9.564397685542594e-01
    #     # PDLP_data.scaled_problem.constraint_matrix[1,2] = 0.0
    #     PDLP_data.scaled_problem.constraint_matrix[2,1] = 5.544352426367563e-02
    #     PDLP_data.scaled_problem.constraint_matrix[2,2] = 5.604930379483980e-01
    #     PDLP_data.scaled_problem.constraint_matrix[3,1] = 2.906560140932720e-02
    #     PDLP_data.scaled_problem.constraint_matrix[3,2] = -5.684713083385977e-01
    #     PDLP_data.scaled_problem.constraint_matrix[4,1] = 2.633791668062104e-03
    #     PDLP_data.scaled_problem.constraint_matrix[4,2] = -5.770712623857750e-01

    #     PDLP_data.scaled_problem.right_hand_side[1] = 9050.578287200928
    #     PDLP_data.scaled_problem.right_hand_side[2] = 2250.650452541391
    #     PDLP_data.scaled_problem.right_hand_side[3] = -1475.098478543439
    #     PDLP_data.scaled_problem.right_hand_side[4] = -2232.497020421451

    #     PDLP_data.constraint_rescaling[1] = 1.000000000000
    #     PDLP_data.constraint_rescaling[2] = 17.250702967684
    #     PDLP_data.constraint_rescaling[3] = 32.906243882067
    #     PDLP_data.constraint_rescaling[4] = 363.141770152987

    #     PDLP_data.scaled_problem.variable_lower_bounds[1,1] = -Inf
    #     PDLP_data.scaled_problem.variable_lower_bounds[1,2] = 3078.934832240113
    #     PDLP_data.scaled_problem.variable_upper_bounds[1,1] = Inf
    #     PDLP_data.scaled_problem.variable_upper_bounds[1,2] = 3079.432769053647

    #     PDLP_data.variable_rescaling[1,1] = 1.045544144940
    #     PDLP_data.variable_rescaling[1,2] = 627.553596380150

    #     PDLP_data.scaled_problem.objective_vector[1,1] = 0.956439768554
    #     PDLP_data.scaled_problem.objective_vector[1,2] = 0.000000000000

    #     PDLP_data.original_problem.objective_vector[1,1] = 1.00000000000
    #     PDLP_data.original_problem.objective_vector[1,2] = 0.00000000000

    #     PDLP_data.original_problem.right_hand_side[1] = 9050.578287200928
    #     PDLP_data.original_problem.right_hand_side[2] = 38825.302440875545
    #     PDLP_data.original_problem.right_hand_side[3] = -48539.950285015781
    #     PDLP_data.original_problem.right_hand_side[4] = -810712.919857113971

    # end
    
    # println("scaled lower variable bounds: \n")
    # println(PDLP_data.scaled_problem.variable_lower_bounds)
    # println("scaled upper variable bounds: \n")
    # println(PDLP_data.scaled_problem.variable_upper_bounds)
    # # Scale the primal weight if desired (otherwise it should be 1.0)
    println("calculating primal weight...")
    if PDLP_data.parameters.bound_objective_rescaling

        PDLP_data.primal_weight .= 1.0

    elseif PDLP_data.parameters.scale_initial_primal_weight

        select_initial_primal_weight(PDLP_data.primal_weight, PDLP_data.original_problem, PDLP_data.dims)
    
    end
    # Come up with a starting step size (Could also put this inside the kernel)
    println("calculating step size...")
    # println("guess vector = $(PDLP_data.kernel_storage.eigenvector)")
    # println("new vector = $(PDLP_data.kernel_storage.new_eigenvector)")
    # println("u vector = $(PDLP_data.kernel_storage.u_vector)")
    update_constant_step_size(PDLP_data.scaled_problem, 
                              PDLP_data.step_size, 
                              PDLP_data.kernel_storage.eigenvector, 
                              PDLP_data.kernel_storage.new_eigenvector, 
                              PDLP_data.kernel_storage.u_vector,
                              PDLP_data.sparsity.nz_count[PDLP_data.dims.current_LP_length],
                              PDLP_data.sparsity.nz_rows,
                              PDLP_data.sparsity.nz_cols,
                              PDLP_data.active_constraint,
                              PDLP_data.dims)
    
    # println("Step size = $(PDLP_data.step_size)")
    # error("avocado!")
    
    # println("Step Size = $(PDLP_data.step_size), Primal Weight = $(PDLP_data.primal_weight)")

    # Run the main loop kernel
    max_size = max(PDLP_data.dims.n_vars, PDLP_data.dims.current_LP_length)
    max_req = Int32(min(256, max(32, ceil(Int, max_size/32)*32))) # number of threads per LP based on LP size

    # Reset total solve and iteration number counters
    PDLP_data.global_counter .= Int32(0)
    PDLP_data.iteration_counter .= Int32(0)
    println("Reached 1: about to get into main loop")    
    
    # Call the main PDLP kernel
    CUDA.@sync @cuda blocks=PDLP_data.dims.n_LPs threads=max_req shmem=max_size*sizeof(Float64) main_loop_kernel(
            solutions,
            objectives,
            PDLP_data.original_problem.variable_lower_bounds,
            PDLP_data.original_problem.variable_upper_bounds,
            PDLP_data.original_problem.right_hand_side,
            PDLP_data.original_problem.objective_vector,
            PDLP_data.original_problem.objective_constant,
            PDLP_data.scaled_problem.variable_lower_bounds,
            PDLP_data.scaled_problem.variable_upper_bounds,
            PDLP_data.scaled_problem.constraint_matrix,
            PDLP_data.scaled_problem.right_hand_side,
            PDLP_data.scaled_problem.objective_vector,
            PDLP_data.sparsity.nz_count[PDLP_data.dims.current_LP_length],
            PDLP_data.sparsity.nz_rows,
            PDLP_data.sparsity.nz_cols,
            PDLP_data.active_constraint,
            PDLP_data.variable_rescaling, 
            PDLP_data.constraint_rescaling, 
            PDLP_data.kernel_storage.current_primal_solution,
            PDLP_data.kernel_storage.current_dual_solution,
            PDLP_data.kernel_storage.current_dual_product,
            PDLP_data.kernel_storage.current_primal_product,
            PDLP_data.kernel_storage.current_primal_gradient,
            PDLP_data.kernel_storage.initial_primal_solution,
            PDLP_data.kernel_storage.initial_dual_solution,
            PDLP_data.kernel_storage.pdhg_primal_solution,
            PDLP_data.kernel_storage.pdhg_dual_solution,
            PDLP_data.kernel_storage.reflected_primal_solution, # 10
            PDLP_data.kernel_storage.reflected_dual_solution, # 11
            PDLP_data.kernel_storage.dual_slack, # 12
            PDLP_data.kernel_storage.residual_primal_product, 
            PDLP_data.kernel_storage.residual_dual_product,
            PDLP_data.kernel_storage.primal_residual,
            PDLP_data.kernel_storage.primal_slack,
            PDLP_data.kernel_storage.dual_residual,
            PDLP_data.kernel_storage.original_primal_solution,
            PDLP_data.kernel_storage.original_primal_gradient,
            PDLP_data.kernel_storage.original_dual_solution,
            PDLP_data.kernel_storage.original_primal_product,
            PDLP_data.kernel_storage.buffer_kkt_primal_solution,
            PDLP_data.kernel_storage.buffer_kkt_primal_product,
            PDLP_data.kernel_storage.buffer_kkt_lower_variable_violation,
            PDLP_data.kernel_storage.buffer_kkt_upper_variable_violation,
            PDLP_data.kernel_storage.buffer_kkt_reduced_costs,
            PDLP_data.kernel_storage.delta_primal,
            PDLP_data.kernel_storage.delta_primal_product,
            PDLP_data.kernel_storage.delta_dual,
            PDLP_data.primal_weight,
            PDLP_data.step_size,
            PDLP_data.termination_reason,
            PDLP_data.dims.current_LP_length,
            PDLP_data.dims.total_LP_length,
            PDLP_data.dims.n_LPs,
            PDLP_data.dims.n_vars,
            PDLP_data.parameters.iteration_limit,
            PDLP_data.parameters.kkt_matrix_pass_limit,
            PDLP_data.parameters.termination_evaluation_frequency,
            PDLP_data.parameters.necessary_reduction_for_restart,
            PDLP_data.parameters.sufficient_reduction_for_restart,
            PDLP_data.parameters.artificial_ratio_for_restart,
            PDLP_data.parameters.extrapolation_coefficient,
            PDLP_data.parameters.reflection_coefficient,
            PDLP_data.parameters.pid_KP,
            PDLP_data.parameters.pid_KI,
            PDLP_data.parameters.pid_KD,
            PDLP_data.parameters.i_smooth,
            PDLP_data.parameters.termination_criteria.eps_optimal_absolute,
            PDLP_data.parameters.termination_criteria.eps_optimal_relative,
            PDLP_data.parameters.termination_criteria.eps_primal_infeasible,
            PDLP_data.parameters.termination_criteria.eps_dual_infeasible,
            return_code,
            global_upper_bound,
            PDLP_data.parameters.skip_hard_problems,
            PDLP_data.global_counter,
            PDLP_data.iteration_counter,
            PDLP_data.skip_flag,
            PDLP_data.iterations,
            PDLP_data.kernel_storage.temp,
            PDLP_data.kernel_storage.temp_dual,
            PDLP_data.kernel_storage.residual_delta_dual,
            )
    return nothing
end

# Check that all the LPs we've created are valid and that the sizes of objects are correct
function validate(PDLP_data::PDLPData)
    # Make sure all entries in the constraint matrix and right-hand side are valid
    if any(!isfinite, PDLP_data.original_problem.constraint_matrix)
        error("Something in the constraint matrix is NaN or Inf")
    end
    if any(!isfinite, PDLP_data.original_problem.right_hand_side)
        error("Something in the right-hand side is NaN or Inf")
    end

    # Make sure the dimensions all line up
    n_LPs = PDLP_data.dims.n_LPs
    n_vars = PDLP_data.dims.n_vars
    tot_len = PDLP_data.dims.total_LP_length

    if size(PDLP_data.original_problem.variable_lower_bounds, 1) < n_LPs || 
       size(PDLP_data.original_problem.variable_lower_bounds, 2) != n_vars
        @show size(PDLP_data.original_problem.variable_lower_bounds)
        @show (n_LPs, n_vars)
        error("Lower bound matrix is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.variable_upper_bounds, 1) < n_LPs || 
       size(PDLP_data.original_problem.variable_upper_bounds, 2) != n_vars
        @show size(PDLP_data.original_problem.variable_upper_bounds)
        @show (n_LPs, n_vars)
        error("Lower bound matrix is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.constraint_matrix, 1) < n_LPs*tot_len || 
       size(PDLP_data.original_problem.constraint_matrix, 2) != n_vars
        @show size(PDLP_data.original_problem.constraint_matrix)
        @show (n_LPs*tot_len, n_vars)
        error("Constraint matrix is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.right_hand_side, 1) < n_LPs*tot_len
        @show size(PDLP_data.original_problem.right_hand_side)
        @show (n_LPs*tot_len)
        error("Right-hand side is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.objective_vector, 1) < n_LPs || 
       size(PDLP_data.original_problem.objective_vector, 2) != n_vars
       @show size(PDLP_data.original_problem.objective_vector)
       @show (n_LPs, n_vars)
        error("Objective vector is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.objective_constant, 1) < n_LPs
        @show size(PDLP_data.original_problem.objective_constant)
        @show n_LPs
        error("Objective constant is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    return nothing
end

function rescale_problem(
    original_problem::LinearProgramSet, 
    scaled_problem::LinearProgramSet, 
    variable_rescaling::CuArray{Float64}, 
    constraint_rescaling::CuArray{Float64},  
    dims::PDLPDims,
    params::PDLPParams
    )

    # Copy problem info to the scaled problem so as to not overwrite the original
    copyto!(scaled_problem, original_problem)

    # Reset rescaling values to be 1.0.
    constraint_rescaling .= 1.0
    variable_rescaling .= 1.0

    # Perform L_inf Ruiz rescaling for `ruiz_iterations` iterations
    if params.ruiz_iterations > 0 # Default is 10
        ruiz_rescaling(
            scaled_problem, 
            params.ruiz_iterations, 
            variable_rescaling, 
            constraint_rescaling, 
            dims,
            )
    end

    # Perform Pock-Chambolle rescaling if alpha isn't nothing
    if !isnothing(params.pock_chambolle_alpha) # Default is 1.0
        pock_chambolle_rescaling(
            scaled_problem, 
            params.pock_chambolle_alpha, 
            variable_rescaling, 
            constraint_rescaling,  
            dims,
            )
    end

    
    return nothing
end

function reset_all_fields!(PDLP_data::PDLPData)

    # Reset all kernel storage
    for field in fieldnames(KernelStorage)
        
        if field == :eigenvector
            CUDA.fill!(getfield(PDLP_data.kernel_storage, field), 1.0)# same vector for each LP
            
        else
            CUDA.fill!(getfield(PDLP_data.kernel_storage, field), 0.0)
        end
        
    end

    # Reset the termination indicator
    CUDA.fill!(PDLP_data.termination_reason, TERMINATION_REASON_UNSPECIFIED)
end

# function calc_problem_norms(problem::LinearProgramSet, dims::PDLPDims)

#     right_hand_side = problem.right_hand_side
#     obj_vec = problem.objective_vector
    
#     const_norm = problem.constraint_bound_norm
#     obj_vec_norm = problem.objective_vector_norm

#     n_vars = dims.n_vars

#     # Identify the number of blocks to use
#     GPU_blocks = Int32(CUDA.attribute(CUDA.device(), CUDA.DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT))

#     CUDA.@sync @cuda blocks=GPU_blocks threads=512 get_norm(right_hand_side, const_norm)

#     CUDA.@sync @cuda blocks=GPU_blocks threads=512 get_norm(obj_vec, obj_vec_norm)

    
#     return nothing
# end

# function get_norm(arr::CuArray{Float64}, norm::Float64)
#     idx = threadIdx().x
# end
include("reduced_load_run_utils.jl")

# The following code runs ucmp for a smaller number of iterations, then solve mpacopf with the binary solution fixed
function solve_ucmp_then_fix_and_solve_mpacopf(file, load_file, gen_prefix)
    use_gpu = false
    env, ucmpmodel = solve_ucmp_reduced_loads(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=100, inner_iterlim=200, start_period=1, end_period=6, scale=1e-4, tight_factor=0.99, use_gpu=use_gpu, verbose=1);
    uc_sol = Int.(round.(ucmpmodel.uc_solution.u_curr))
    ntime = div(size(uc_sol)[2], 3)
    on_status = Vector{Int}[uc_sol[:,3*i-2] for i in 1:ntime]
    switch_on = Vector{Int}[uc_sol[:,3*i-1] for i in 1:ntime]
    switch_off = Vector{Int}[uc_sol[:,3*i] for i in 1:ntime]
    total_con = sum(sum(ucmpmodel.uc_params.con .* switch_on[i]) for i in 1:ntime)
    total_coff = sum(sum(ucmpmodel.uc_params.coff .* switch_off[i]) for i in 1:ntime)
    env, mpmodel = solve_mpacopf_reduced_loads(file, load_file; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, outer_iterlim=1000, inner_iterlim=200, start_period=1, end_period=6, scale=1e-4, tight_factor=0.99, use_gpu=use_gpu, on_status=on_status, verbose=1);
    
    println("Total on cost: $(total_con)")
    println("Total off cost: $(total_coff)")
    println("UC cost + mpmodel cost: $(total_con + total_coff + mpmodel.info.objval)")
    println("ucmpmodel objective value: $(ucmpmodel.info.objval)")
end

case = "case118"
file = "/home/wzhang483/gpu_acopf/data/$(case).m"
load_file = "/home/wzhang483/gpu_acopf/data/multiperiod_data/$(case)_onehour_168"
gen_prefix = "/home/wzhang483/gpu_acopf/data/multiperiod_data/$(case)_gen"

solve_ucmp_then_fix_and_solve_mpacopf(file, load_file, gen_prefix)

# To directly solve ucmp, run as follows with large enough outer_iterlim and inner_iterlim
# env, ucmpmodel = solve_ucmp_reduced_loads(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=1000, inner_iterlim=200, start_period=1, end_period=6, scale=1e-4, tight_factor=0.99, use_gpu=use_gpu, verbose=1);

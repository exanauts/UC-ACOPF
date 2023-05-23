include("reduced_load_run_utils.jl")

case = ARGS[1]
T = parse(Int, ARGS[2])
file = "./data/$(case).m"
load_file = "./data/multiperiod_data/$(case)_onehour_168"
gen_prefix = "./data/multiperiod_data/$(case)_gen_initial_all_on"

use_gpu = false
if length(ARGS) > 2 && ARGS[3] == "use_gpu"
    use_gpu = true
end
env, ucmpmodel = solve_ucmp_reduced_loads(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e4, rho_uc=1e4, outer_iterlim=50, inner_iterlim=200, start_period=1, end_period=T, scale=1e-4, tight_factor=0.99, use_gpu=use_gpu, verbose=1);

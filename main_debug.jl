using ExaAdmm
using CUDA

case = "case9"
# case = "case118"
# case = "case2868rte"
file = "./data/$(case).m"
load_file = "./data/multiperiod_data/$(case)_onehour_168"
gen_prefix = "./data/multiperiod_data/$(case)_gen"

# env, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.03, rho_pq=1e3, rho_va=1e4, rho_uc=1e4, outer_iterlim=6, inner_iterlim=50, start_period=1, end_period=2, scale=1e-4, tight_factor=0.99, use_gpu=false, warm_start=false); # FOR GPU DEBUGGING
# env, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.03, rho_pq=1e3, rho_va=1e4, rho_uc=1e4, outer_iterlim=6, inner_iterlim=50, start_period=1, end_period=2, scale=1e-4, tight_factor=0.99, use_gpu=true, warm_start=false); # FOR GPU DEBUGGING

_, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=1, inner_iterlim=5, start_period=1, end_period=6, scale=1e-4, tight_factor=0.99, use_gpu=false, warm_start=false);
_, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=1, inner_iterlim=5, start_period=1, end_period=6, scale=1e-4, tight_factor=0.99, use_gpu=true, warm_start=false);

# Case 2868 rte, 3 hours
# _, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=1, inner_iterlim=5, start_period=1, end_period=3, scale=1e-4, tight_factor=0.99, use_gpu=false, warm_start=false);
# _, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=1, inner_iterlim=5, start_period=1, end_period=3, scale=1e-4, tight_factor=0.99, use_gpu=true, warm_start=false);

println("Run ended.")
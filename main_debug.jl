using ExaAdmm
using CUDA

include("jump_models.jl")

case = "case9"
# case = "case118"
# case = "case2868rte"
file = "./data_v2/$(case).m"
load_file = "./data_v2/multiperiod_data/$(case)_onehour_168"
gen_prefix = "./data_v2/multiperiod_data/$(case)_gen"

T = 24

# env, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.03, rho_pq=1e3, rho_va=1e4, rho_uc=1e4, outer_iterlim=1000, inner_iterlim=200, start_period=1, end_period=2, scale=1e-4, tight_factor=0.99, use_gpu=false, warm_start=false); # FOR GPU DEBUGGING
# env, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.03, rho_pq=1e3, rho_va=1e4, rho_uc=1e4, outer_iterlim=1000, inner_iterlim=200, start_period=1, end_period=2, scale=1e-4, tight_factor=0.99, use_gpu=true, warm_start=false); # FOR GPU DEBUGGING

# _, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=1, inner_iterlim=5, start_period=1, end_period=6, scale=1e-4, tight_factor=0.99, use_gpu=false, warm_start=false);
# _, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=1, inner_iterlim=5, start_period=1, end_period=6, scale=1e-4, tight_factor=0.99, use_gpu=true, warm_start=false);

# env, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e5, outer_iterlim=1000, inner_iterlim=200, start_period=1, end_period=6, scale=1e-4, tight_factor=0.99, use_gpu=true, warm_start=true);

# Case 2868 rte, 3 hours
# _, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=1, inner_iterlim=5, start_period=1, end_period=3, scale=1e-4, tight_factor=0.99, use_gpu=false, warm_start=false);
# _, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e5, rho_uc=1e3, outer_iterlim=10, inner_iterlim=100, start_period=1, end_period=3, scale=1e-4, tight_factor=0.99, use_gpu=true, warm_start=false);

env, ucmpmodel = ExaAdmm.solve_ucmp(file, load_file, gen_prefix; ramp_ratio=0.1, rho_pq=1e3, rho_va=1e4, rho_uc=1e4, outer_iterlim=100, inner_iterlim=50, start_period=1, end_period=T, scale=1e-4, tight_factor=0.99, use_gpu=false, warm_start=true);

u_on = zeros(ucmpmodel.mpmodel.models[1].grid_data.ngen, T)
for g in 1:ucmpmodel.mpmodel.models[1].grid_data.ngen
    for t in 1:T
        u_on[g,t] = ucmpmodel.uc_solution.u_curr[g,3*t-2]
    end
end
# solve_multiperiod_with_uc(circuit, load, T, u_on, u_su, u_sd)


println("Run ended.")
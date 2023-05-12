using ExaAdmm
using CUDA

function solve_mpacopf_reduced_loads(case::String, load_prefix::String;
    case_format="matpower", start_period=1, end_period=1,
    outer_iterlim=20, inner_iterlim=1000, rho_pq=4e2, rho_va=4e4,
    obj_scale=1.0, scale=1e-4,
    use_gpu=false, use_linelimit=true, tight_factor=1.0,
    outer_eps=2*1e-4, gpu_no=0, verbose=1, ramp_ratio=0.02, warm_start=true, multiperiod_tight=true,
    on_status=nothing, switch_on=nothing, switch_off=nothing)

    T = Float64; TD = Array{Float64,1}; TI = Array{Int,1}; TM = Array{Float64,2}
    if use_gpu
        CUDA.device!(gpu_no)
        TD = CuArray{Float64,1}; TI = CuArray{Int,1}; TM = CuArray{Float64,2}
    end

    env = ExaAdmm.AdmmEnv{T,TD,TI,TM}(case, rho_pq, rho_va; case_format=case_format,
            use_gpu=use_gpu, use_linelimit=use_linelimit,
            load_prefix=load_prefix, tight_factor=tight_factor, gpu_no=gpu_no, verbose=verbose)

    env.load.pd ./= 8
    env.load.qd ./= 8

    mod = ExaAdmm.ModelMpacopf{T,TD,TI,TM}(env; start_period=start_period, end_period=end_period, ramp_ratio=ramp_ratio, on_status=on_status, switch_on=switch_on, switch_off=switch_off)

    n = mod.models[1].n
    env.params.scale = scale
    env.params.obj_scale = obj_scale
    env.params.outer_eps = outer_eps
    env.params.outer_iterlim = outer_iterlim
    env.params.inner_iterlim = inner_iterlim

    if warm_start
        env.params.verbose = 0
        for i=1:mod.len_horizon
            ExaAdmm.admm_two_level(env, mod.models[i])
        end
        env.params.verbose = 1
        ExaAdmm.init_solution!(mod, mod.solution, rho_pq, rho_va)
    end

    ExaAdmm.admm_two_level(env, mod)

    return env, mod
end

function solve_ucmp_reduced_loads(case::String, load_prefix::String, gen_prefix::String;
    case_format="matpower", start_period=1, end_period=1,
    outer_iterlim=20, inner_iterlim=1000, rho_pq=400.0, rho_va=40000.0, rho_uc=40000.0,
    obj_scale=1.0, scale=1e-4, storage_ratio=0.0, storage_charge_max=1.0,
    use_gpu=false, ka_device=nothing, use_linelimit=true, use_projection=false, tight_factor=1.0,
    outer_eps=2*1e-4, gpu_no=0, verbose=1, ramp_ratio=0.02, warm_start=true, multiperiod_tight=true
)

    T = Float64
    TD = Array{Float64,1}; TI = Array{Int,1}; TM = Array{Float64,2}
    ka_device = nothing

    
    env = ExaAdmm.AdmmEnv{T,TD,TI,TM}(case, rho_pq, rho_va; load_prefix=load_prefix, case_format=case_format,
            use_gpu=use_gpu, ka_device=ka_device, use_linelimit=use_linelimit,
            use_projection=use_projection, tight_factor=tight_factor, gpu_no=gpu_no,
            storage_ratio=storage_ratio, storage_charge_max=storage_charge_max,
            verbose=verbose)
            
    env.load.pd ./= 8
    env.load.qd ./= 8
        
    mod = ExaAdmm.UCMPModel{T,TD,TI,TM}(env, gen_prefix; start_period=start_period, end_period=end_period, ramp_ratio=ramp_ratio)

    env.params.scale = scale
    env.params.obj_scale = obj_scale
    env.params.outer_eps = outer_eps
    env.params.outer_iterlim = outer_iterlim
    env.params.inner_iterlim = inner_iterlim

    if warm_start
        env.params.verbose = 0
        ExaAdmm.admm_two_level(env, mod.mpmodel)
        env.params.verbose = 1
        ExaAdmm.init_solution!(mod, mod.uc_solution, rho_pq, rho_va, rho_uc)
    end

    ExaAdmm.admm_two_level(env, mod)

    return env, mod
end
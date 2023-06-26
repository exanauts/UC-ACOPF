using Ipopt
using Juniper
using Ipopt, HiGHS

include("MPOPF-sampler/src/mpc_data.jl")
include("jump_models.jl")

# case = "case9"
case = "case118"
# case = ARGS[1]
gen_scenario = "$(case)_gen"

T = 24
# T = parse(Int, ARGS[2])

println("============================================")
println("Solving $(case) with $(T) hours...")

circuit = getcircuit("data_v2/$(case)", 100, 0.1)
load = getload("data_v2/multiperiod_data/$(case)_onehour_168", circuit)

v0 = readdlm("data_v2/multiperiod_data/$(gen_scenario).v0", Int)
hu = readdlm("data_v2/multiperiod_data/$(gen_scenario).hu", Int)
hd = readdlm("data_v2/multiperiod_data/$(gen_scenario).hd", Int)
tu = readdlm("data_v2/multiperiod_data/$(gen_scenario).tu", Int)
td = readdlm("data_v2/multiperiod_data/$(gen_scenario).td", Int)
con = readdlm("data_v2/multiperiod_data/$(gen_scenario).con")
coff = readdlm("data_v2/multiperiod_data/$(gen_scenario).coff")


model = get_ucmodel(circuit, load, T, v0, tu, td, hu, hd, con, coff)

nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
# minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver)
mip_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver, "mip_solver"=>mip_solver)


set_optimizer(model, minlp_solver)

@time optimize!(model)

println("Termination status: ", termination_status(model))
println("============================================")
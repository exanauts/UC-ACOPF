using Ipopt
using Juniper
using Ipopt, HiGHS

include("MPOPF-sampler/src/mpc_data.jl")
include("jump_models.jl")

case = "case9"
# case = "case118"
gen_scenario = "$(case)_gen"

T = 6
circuit = getcircuit("data/$(case)", 100, 0.1)
load = getload("data/multiperiod_data/$(case)_onehour_168", circuit)
# load.pd ./= 3
# load.qd ./= 3

v0 = readdlm("data/multiperiod_data/$(gen_scenario).v0", Int)
hu = readdlm("data/multiperiod_data/$(gen_scenario).hu", Int)
hd = readdlm("data/multiperiod_data/$(gen_scenario).hd", Int)
tu = readdlm("data/multiperiod_data/$(gen_scenario).tu", Int)
td = readdlm("data/multiperiod_data/$(gen_scenario).td", Int)
con = readdlm("data/multiperiod_data/$(gen_scenario).con")
coff = readdlm("data/multiperiod_data/$(gen_scenario).coff")


model = get_ucmodel(circuit, load, T, v0, tu, td, hu, hd, con, coff)

nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
# minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver)
mip_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver, "mip_solver"=>mip_solver)


set_optimizer(model, minlp_solver)

optimize!(model)
using JuMP, Ipopt

include("MPOPF-sampler/src/mpc_data.jl")
include("jump_models.jl")

case = "case118"
T = 24

println("============================================")
println("Solving continuous relaxation with Ipopt for $(case) with $(T) hours...")

circuit = getcircuit("data_v2/$(case)", 100, 0.1)
load = getload("data_v2/multiperiod_data/$(case)_onehour_168", circuit)

gen_scenario = "$(case)_gen"
v0 = readdlm("data_v2/multiperiod_data/$(gen_scenario).v0", Int)
hu = readdlm("data_v2/multiperiod_data/$(gen_scenario).hu", Int)
hd = readdlm("data_v2/multiperiod_data/$(gen_scenario).hd", Int)
tu = readdlm("data_v2/multiperiod_data/$(gen_scenario).tu", Int)
td = readdlm("data_v2/multiperiod_data/$(gen_scenario).td", Int)
con = readdlm("data_v2/multiperiod_data/$(gen_scenario).con")
coff = readdlm("data_v2/multiperiod_data/$(gen_scenario).coff")

model = get_ucmodel(circuit, load, T, v0, tu, td, hu, hd, con, coff, cont_relax=true)
optimizer = optimizer_with_attributes(
    Ipopt.Optimizer, 
    "linear_solver" => "ma27",
    "print_level" => 0
)
set_optimizer(model, optimizer)
time = @elapsed optimize!(model)

println("Termination status: $(termination_status(model))")
println("Objective value: $(objective_value(model))")
println("Solution time: $(time)")
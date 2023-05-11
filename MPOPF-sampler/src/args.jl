using DataStructures

#                  ramp scale (original demand)  After random demand
# case9              0.0004                        0.03
# case30             0.0005                        0.03
# case118            0.002                         0.04
# case300            0.002                         0.06
# case1354pegase     0.0025
# case2383wp         0.005
# case9241pegase     0.01
case9_args = OrderedDict{String, Any}(
    "case" => "data/case9",
    "scen" => "data/case9/halfhour_30",
    "T" => 30,
    "H" => 1,
    "load_scale" => 0.1,
    "ramp_scale" => 0.03,
    "warm" => "warm",
    "opt" => 2,
    "cut_line" => true,
    "cut_gen" => false,
    "perturb" => 0,
    "qp" => 0,
    "load_sol" => 0,
    "powerflow_solve" => false,
    "profname" => nothing,
    "cut_time" => 1
)

case30_args = OrderedDict{String, Any}(
    "case" => "data/case30",
    "scen" => "data/case30/halfhour_30",
    "T" => 30,
    "H" => 1,
    "load_scale" => 0.05,
    "ramp_scale" => 0.04,
    "warm" => "warm",  # FIXME This doesn't work well
#    "warm" => "cold",
    "opt" => 2,
    "cut_line" => true,
#    "cut_line" => false,
    "cut_gen" => false,
    "perturb" => 0,
    "qp" => 0,
    "load_sol" => 0,
    "powerflow_solve" => false,
    "profname" => nothing,
    "cut_time" => 1
)

case118_args = OrderedDict{String, Any}(
    "case" => "data/case118",
    "scen" => "data/case118/halfhour_30",
    "T" => 30,
    "H" => 1,
    "load_scale" => 0.1,
    "ramp_scale" => 0.04,
    "warm" => "warm",
    "opt" => 2,
    "cut_line" => true,
    "cut_gen" => false,
    "perturb" => 0,
    "qp" => 0,
    "load_sol" => 0,
    "powerflow_solve" => false,
    "profname" => nothing,
    "cut_time" => 1
)

case300_args = OrderedDict{String, Any}(
    "case" => "data/case300",
    "scen" => "data/case300/halfhour_30",
    "T" => 30,
    "H" => 1,
    "load_scale" => 0.1,
    "ramp_scale" => 0.06,
    "warm" => "cold",
    "opt" => 2,
    "cut_line" => true,
    "cut_gen" => false,
    "perturb" => 0,
    "qp" => 0,
    "load_sol" => 0,
    "powerflow_solve" => false,
    "profname" => nothing,
    "cut_time" => 1
)

case1354pegase_args = OrderedDict{String, Any}(
    "case" => "data/case1354pegase",
    "scen" => "data/case1354pegase/halfhour_30",
    "T" => 10,
    "H" => 1,
    "load_scale" => 0.1,
    "ramp_scale" => 0.08,
    "warm" => "warm",
#    "warm" => "cold",
    "opt" => 2,
    "cut_line" => false,
    "cut_gen" => false,
    "perturb" => 0,
    "qp" => 0,
    "load_sol" => 0,
    "powerflow_solve" => false,
    "profname" => nothing,
    "cut_time" => 1
)


#!/bin/zsh
#!/bin/bash

# For git bash
#$HOME/AppData/Local/Julia-1.3.1/bin/julia --project=. src/mpc.jl "data/case9" "data/case9/halfhour_30" 10 1 1.0 0.01 "cold" 2 0 1 0 0 0 0 "case9_result"
#                  ramp scale
# case9              0.0004
# case30             0.0005
# case118            0.002
# case300            0.002
# case1354pegase     0.0025
# case2383wp         0.005
# case9241pegase     0.01
#
# Args:
#
#    case = args[1]
#    scen = args[2]
#    T = max(parse(Int,args[3]),1)
#    H = max(parse(Int,args[4]),1)
#    load_scale = parse(Float64,args[5])
#    ramp_scale = parse(Float64,args[6])
#    warm = args[7]
#    opt = max(parse(Int,args[8]),2)
#    cut_line = (parse(Int,args[9]) == 1) ? true : false
#    cut_gen = (parse(Int,args[10]) == 1) ? true : false
#    perturb = parse(Float64,args[11])
#    qp = parse(Int,args[12])
#    load_sol = parse(Int,args[13])
#    powerflow_solve = (parse(Int,args[14]) == 1) ? true : false
#    profname = nothing
#    Optional arg: niter => run niter times on this load_scale setting
#
#    saved dir = "./results/$(basename(case))/"
#
#
#
# For ubuntu zsh
# For cut line || cut gen
#/bin/julia --project=. src/mpc.jl "data/case9" "data/case9/halfhour_30" 10 1 1.0 0.01 "cold" 2 0 1 0 0 0 0 "case9_result"
#
# No cut
#julia --project=. src/mpc.jl "data/case9" "data/case9/halfhour_30" 10 1 0.05 0.0004 "cold" 2 0 0 0 0 0 0 20  > case9_log.txt
julia --project=. src/mpc.jl "data/case9" "data/case9/halfhour_30" 10 1 0.05 0.004 "cold" 2 0 0 0 0 0 0 1  > results/case9_log.txt
#/bin/julia --project=. src/mpc.jl "data/case30" "data/case30/halfhour_30" 30 1 0 0.0005 "cold" 2 0 0 0 0 0 0 10 >> case30_log.txt
#/bin/julia --project=. src/mpc.jl "data/case118" "data/case118/halfhour_30" 30 1 0 0.002 "cold" 2 0 0 0 0 0 0 10 >> case118_log.txt
#/bin/julia --project=. src/mpc.jl "data/case300" "data/case300/halfhour_30" 30 1 0 0.002 "cold" 2 0 0 0 0 0 0 2000
#/bin/julia --project=. src/mpc.jl "data/case1354pegase" "data/case1354pegase/halfhour_30" 10 1 0 0.0025 "cold" 2 0 0 0 0 0 0 2000
#/bin/julia --project=. src/mpc.jl "data/case2383wp" "data/case2383wp/halfhour_30" 10 1 0 0.005 "cold" 2 0 0 0 0 0 0 2000


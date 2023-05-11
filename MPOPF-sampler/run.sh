#!/bin/zsh

#                  ramp scale
# case9              0.0004
# case30             0.0005
# case118            0.002
# case300            0.002
# case1354pegase     0.0025
# case2383wp         0.005
# case9241pegase     0.01
#
#
#

#julia --project=. src/run_mpc.jl "case9" 10000 0.1 0.03 > results/case9_grid.txt
#$HOME/downloads/julia-1.5.0/bin/julia --project=. src/run_mpc.jl "case9" 100 -e 'using Pkg; Pkg.instantiate()' &
#$HOME/downloads/julia-1.5.0/bin/julia --project=. src/run_mpc.jl "case30" 100 -e 'using Pkg; Pkg.instantiate()' &
#$HOME/downloads/julia-1.5.0/bin/julia --project=. src/run_mpc.jl "case118" 100 -e 'using Pkg; Pkg.instantiate()' &
$HOME/downloads/julia-1.5.0/bin/julia --project=. src/run_mpc.jl "case300" 100 -e 'using Pkg; Pkg.instantiate()'
#$HOME/downloads/julia-1.5.0/bin/julia --project=. src/run_mpc.jl "case1354pegase" 100 -e 'using Pkg; Pkg.instantiate()'


#                                         10 -> number of MPOPF problems that are solved (different load demands)
#julia --project=. src/run_mpc.jl "case9" 10 -e 'using Pkg; Pkg.instantiate()'

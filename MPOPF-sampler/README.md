# MPOPF Sampler

This is an implementation of a multi-period AC optimal power flow (MPOPF) problem with N-1 line contingency and random load demands.

# How to run (need julia to be called by command `julia`)

```shell
$ sh run.sh
```

After running this code, results will be saved in the `(PROJECT)/results` folder.
In the `info_ipopt_caseXXXX.jld2` file, values of $(P_g, Q_g, V_m, V_a, P_d, Q_d)$ for the whole period are saved.
Sample extraction code is given in `extract.jl`. Run `julia extract.jl` for an example.

In the following code from `run.sh`, you may substitute `case9` with `case9, case30, case118, case300`.
The next argument (Integer) is the number of MPOPF problems to be solved under different load demands.

```
julia --project=. src/run_mpc.jl "case9" 10 -e 'using Pkg; Pkg.instantiate()'
```

To change other options, you can change it from `src/args.jl` where

* case: str, the name of the case file
* scen: str, the name of the scenario file (load profile)
* T: int, the length of a time horizon, e.g., T=10 for 10 time periods in a horizon.
* load_scale: float, load scale in (0.0,1.0], e.g., 0.1 means 10% uniform perturbation in load
* ramp_scale: float, ramp scale in (0.0,1.0], e.g., RS=0.01 means that we allow generator ramp rate to be a 1 percent.
* warm: str in {"cold", "shift_copy", "shift_phase1"}. This is for warm-starting of a rolling horizon.
* opt: int, the option file number IPOPT will be using. OPT=2 will use the default file ipopt.opt.
* cut_line: bool, true if we want to cut a line, false otherwise.
* cut_time: time when line contingency happens

If you want to change your model for line contingency, you can do so by changing the function in line 1084 `get_cut_mpcmodel` which are models from `mpc_cut_model.jl`.

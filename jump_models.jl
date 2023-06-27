using JuMP
using DelimitedFiles
function get_mpcmodel(circuit, demand, T; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing, u_on=nothing, u_su=nothing, u_sd=nothing, con=nothing, coff=nothing)
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    m = Model()

    # Shortcuts
    baseMVA = circuit.baseMVA
    busref = circuit.busref
    bus = circuit.bus
    line = circuit.line
    gen = circuit.gen
    yline = circuit.yline
    ybus = circuit.ybus
    busdict = circuit.busdict
    frombus = circuit.frombus
    tobus = circuit.tobus
    bus2gen = circuit.bus2gen

    Pd = demand.pd
    Qd = demand.qd
    # T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)
    num_lines = length(line)

    if isnothing(u_on)
        u_on = ones(Int, num_gens, T)
        u_su = zeros(Int, num_gens, T)
        u_sd = zeros(Int, num_gens, T)
    end
    

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    # @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    # @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, gen[g].Pmin * u_on[g,t] <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax * u_on[g,t])
    @variable(m, gen[g].Qmin * u_on[g,t] <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax * u_on[g,t])
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        set_lower_bound(Va[t,busref], bus[busref].Va)
        set_upper_bound(Va[t,busref], bus[busref].Va)
    end

    #= Objective function
    If not piecewise, then it is a quadratic sum of the real power at (t,g)
    \sum_t\sum_g gencost1[g] * (baseMVA * Pg[t,g])^2
            + gencost2[g] * (baseMVA * Pg[t,g])
            + gencost3[g]
    =#
    uc_cost = sum(con[g] * sum(u_su[g,t] for t=1:T) + coff[g] * sum(u_sd[g,t] for t=1:T) for g=1:num_gens)
    if piecewise
        @variable(m, Cg[t=1:T,g=1:num_gens])
        @NLobjective(m, Min, sum(Cg[t,g] for t=1:T,g=1:num_gens) + uc_cost)
        @constraint(m, plcurve[t=1:T,g=1:num_gens,p=1:gen[g].n-1],
		    Cg[t,g] - (((gen[g].coeff[2*p+2] - gen[g].coeff[2*p])/(gen[g].coeff[2*p+1] - gen[g].coeff[2*p-1]))*(baseMVA*Pg[t,g] - gen[g].coeff[2*p-1]) + gen[g].coeff[2*p]) >= 0
                   )
    else
        @NLobjective(m, Min, sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2
                                + gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g])
                                + gen[g].coeff[gen[g].n] for t=1:T,g=1:num_gens) + uc_cost)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            # @constraint(m, ramping[t=1:T-1,g=1:num_gens],
            #         -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
            @constraint(m, ramping[t=1:T-1,g=1:num_gens],
                    -gen[g].ramp_agc * u_on[g,t+1] - gen[g].Pmax * u_sd[g,t+1] <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc * u_on[g,t] + gen[g].Pmax * u_su[g,t+1])
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    # Power flow constraints: real part
    @NLconstraint(m, pfreal[t=1:T,b=1:num_buses],
                  (sum(yline[l].YffR for l in frombus[b])
                   + sum(yline[l].YttR for l in tobus[b])
                   + ybus[b].YshR)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfI*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA
                  == 0)

    # Power flow constraints: imaginary part
    @NLconstraint(m, pfimag[t=1:T,b=1:num_buses],
                  (sum(-yline[l].YffI for l in frombus[b])
                   + sum(-yline[l].YttI for l in tobus[b])
                   - ybus[b].YshI)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (-yline[l].YftI*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftR*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (-yline[l].YtfI*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfR*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA
                  == 0)

    # Line limits
    # rateA - admittance rate?
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Yre = zeros(num_linelimits)
    Yim = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2
        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR
    end

    @NLconstraint(m, flowmaxfrom[t=1:T,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].from]]^2 *
                  (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    - Yim[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <= 0)

    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (to bus)
        l = limind[i]
        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end

    @NLconstraint(m, flowmaxto[t=1:T,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].to]]^2 *
                  (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    -Yim[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <=0)
    return m
end

function get_ucmodel(circuit, demand, T,
                     v0, tu, td, hu, hd, con, coff; 
                     has_ramping=true,
                     phase1=false, piecewise=false,
                     prev_val=nothing, cont_relax=false)
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    m = Model()

    # Shortcuts
    baseMVA = circuit.baseMVA
    busref = circuit.busref
    bus = circuit.bus
    line = circuit.line
    gen = circuit.gen
    yline = circuit.yline
    ybus = circuit.ybus
    busdict = circuit.busdict
    frombus = circuit.frombus
    tobus = circuit.tobus
    bus2gen = circuit.bus2gen

    Pd = demand.pd
    Qd = demand.qd
    # T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)
    num_lines = length(line)
    
    #=
    u_on[g,t]: generator g is on or off at time t
    u_su[g,t]: generator g is started up at time t
    u_sd[g,t]: generator g is shut down at time t
    =#
    if cont_relax
        @variable(m, u_on[g=1:num_gens,t=1:T], lower_bound=0, upper_bound=1)
        @variable(m, u_su[g=1:num_gens,t=1:T], lower_bound=0, upper_bound=1)
        @variable(m, u_sd[g=1:num_gens,t=1:T], lower_bound=0, upper_bound=1)
    else
        @variable(m, u_on[g=1:num_gens,t=1:T], Bin)
        @variable(m, u_su[g=1:num_gens,t=1:T], Bin)
        @variable(m, u_sd[g=1:num_gens,t=1:T], Bin)
    end
    @constraint(m, uc_state_init[g=1:num_gens], v0[g] - u_on[g,1] + u_su[g,1] - u_sd[g,1] == 0)
    @constraint(m, uc_state[g=1:num_gens,t=1:T-1], u_on[g,t] - u_on[g,t+1] + u_su[g,t+1] - u_sd[g,t+1] == 0)
    @constraint(m, initial_on[g=1:num_gens], sum(1-u_on[g,t] for t in 1:hu[g]) == 0)
    @constraint(m, initial_off[g=1:num_gens], sum(u_on[g,t] for t in 1:hd[g]) == 0)
    # @constraint(m, min_ontime[g=1:num_gens,t=tu[g]:T], sum(u_su[g,t] for i in t-tu[g]+1:t) <= u_on[g,t])
    # @constraint(m, min_offtime[g=1:num_gens,t=td[g]:T], sum(u_sd[g,t] for i in t-td[g]+1:t) <= 1 - u_on[g,t])
    for g=1:num_gens
        if tu[g] > 0
            @constraint(m, [t=tu[g]:T], sum(u_su[g,t] for i in t-tu[g]+1:t) <= u_on[g,t])
        end
        if td[g] > 0
            @constraint(m, [t=td[g]:T], sum(u_sd[g,t] for i in t-td[g]+1:t) <= 1 - u_on[g,t])
        end
    end
    @expression(m, uc_cost, sum(con[g] * sum(u_su[g,t] for t=1:T) + coff[g] * sum(u_sd[g,t] for t=1:T) for g=1:num_gens) )

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    @variable(m, Pg[t=1:T,g=1:num_gens])
    @variable(m, Qg[t=1:T,g=1:num_gens])
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])

    @constraint(m, pg_lower_bound[t=1:T,g=1:num_gens], Pg[t,g] >= gen[g].Pmin * u_on[g,t])
    @constraint(m, pg_upper_bound[t=1:T,g=1:num_gens], Pg[t,g] <= gen[g].Pmax * u_on[g,t])
    @constraint(m, qg_lower_bound[t=1:T,g=1:num_gens], Qg[t,g] >= gen[g].Qmin * u_on[g,t])
    @constraint(m, qg_upper_bound[t=1:T,g=1:num_gens], Qg[t,g] <= gen[g].Qmax * u_on[g,t])

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        set_lower_bound(Va[t,busref], bus[busref].Va)
        set_upper_bound(Va[t,busref], bus[busref].Va)
    end

    #= Objective function
    If not piecewise, then it is a quadratic sum of the real power at (t,g)
    \sum_t\sum_g gencost1[g] * (baseMVA * Pg[t,g])^2
            + gencost2[g] * (baseMVA * Pg[t,g])
            + gencost3[g]
    =#
    if piecewise
        @variable(m, Cg[t=1:T,g=1:num_gens])
        @NLobjective(m, Min, sum(Cg[t,g] for t=1:T,g=1:num_gens) + uc_cost)
        @constraint(m, plcurve[t=1:T,g=1:num_gens,p=1:gen[g].n-1],
		    Cg[t,g] - (((gen[g].coeff[2*p+2] - gen[g].coeff[2*p])/(gen[g].coeff[2*p+1] - gen[g].coeff[2*p-1]))*(baseMVA*Pg[t,g] - gen[g].coeff[2*p-1]) + gen[g].coeff[2*p]) >= 0
                   )
    else
        @NLobjective(m, Min, sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2
                                + gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g])
                                + gen[g].coeff[gen[g].n] for t=1:T,g=1:num_gens) + uc_cost)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            @constraint(m, ramping_up[t=1:T-1,g=1:num_gens],
                    Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc * u_on[g,t] + gen[g].Pmax * u_su[g,t+1])
            @constraint(m, ramping_down[t=1:T-1,g=1:num_gens],
                    Pg[t+1,g] - Pg[t,g] >= -gen[g].ramp_agc * u_on[g,t+1] - gen[g].Pmax * u_sd[g,t+1])
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    # Power flow constraints: real part
    @NLconstraint(m, pfreal[t=1:T,b=1:num_buses],
                  (sum(yline[l].YffR for l in frombus[b])
                   + sum(yline[l].YttR for l in tobus[b])
                   + ybus[b].YshR)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfI*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA
                  == 0)

    # Power flow constraints: imaginary part
    @NLconstraint(m, pfimag[t=1:T,b=1:num_buses],
                  (sum(-yline[l].YffI for l in frombus[b])
                   + sum(-yline[l].YttI for l in tobus[b])
                   - ybus[b].YshI)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (-yline[l].YftI*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftR*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (-yline[l].YtfI*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfR*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA
                  == 0)

    # Line limits
    # rateA - admittance rate?
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Yre = zeros(num_linelimits)
    Yim = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2
        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR
    end

    @NLconstraint(m, flowmaxfrom[t=1:T,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].from]]^2 *
                  (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    - Yim[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <= 0)

    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (to bus)
        l = limind[i]
        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end

    @NLconstraint(m, flowmaxto[t=1:T,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].to]]^2 *
                  (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    -Yim[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <=0)
    return m
end

function solve_multiperiod(circuit, load, T)
    model = get_mpcmodel(circuit, load, T)
    set_optimizer(model, Ipopt.Optimizer)
    optimize!(model)
    return model
end

function solve_multiperiod_with_uc(circuit, load, T, u_on, u_su, u_sd, con, coff)
    model = get_mpcmodel(circuit, load, T, u_on=u_on, u_su=u_su, u_sd=u_sd, con=con, coff=coff)
    set_optimizer(model, Ipopt.Optimizer)
    optimize!(model)
    return model
end
#=
function feasibility_check(v0, tu, td, hu, hd, u_on)
    T = length(u_on)
    u_su = zeros(Int, T)
    u_sd = zeros(Int, T)
    u_su[1] = max(u_on[1] - v0, 0)
    u_sd[1] = max(v0 - u_on[1], 0)
    for t in 2:T
        u_su[t] = max(u_on[t] - u_on[t-1], 0)
        u_sd[t] = max(u_on[t-1] - u_on[t], 0)
    end

    # This check is not necessary, only for debugging purposes
    @assert v0 - u_on[1] + u_su[1] - u_sd[1] == 0
    for t in 2:T
        @assert u_on[t-1] - u_on[t] + u_su[t] - u_sd[t] == 0
    end

    if sum(1 .- u_on[1:hu]) != 0
        return 0
    end

    if sum(u_on[1:hd]) != 0
        return 0
    end

    for t in tu:T
        if sum(u_su[t-tu+1:t]) > u_on[t]
            return 0
        end
    end

    for t in td:T
        if sum(u_sd[t-td+1:t]) > 1 - u_on[t]
            return 0
        end
    end

    return 1
end

function generate_uc_feasible_solutions(v0, tu, td, hu, hd, T)
    s = Set()
    for i in 0:2^T-1
        u_on = reverse(digits(i, base=2, pad=T))
        if feasibility_check(v0, tu, td, hu, hd, u_on) == 1
            push!(s, i)
        end
    end
    return s
end
=#
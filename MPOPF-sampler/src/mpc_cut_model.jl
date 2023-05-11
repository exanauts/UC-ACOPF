using JuMP
using Ipopt

function get_cut_mpcmodel1(circuit, cut_circuit, demand, cut_time; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing)
    # Model1 Soft constraints all time
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    # TODO: Need to receive cut circuit and cut time and divide it with before cut and after cut
    lambda_p = 3000
    lambda_q = 3000
    lambda_l = 3000

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
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)

    # Line limits
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)
    Yre_from = zeros(num_linelimits)
    Yim_from = zeros(num_linelimits)
    Yre_to = zeros(num_linelimits)
    Yim_to = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2

        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre_from[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim_from[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR

        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre_to[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim_to[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end


    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)

    Yff_abs2_cut = zeros(num_linelimits_cut)
    Yft_abs2_cut = zeros(num_linelimits_cut)
    Ytf_abs2_cut = zeros(num_linelimits_cut)
    Ytt_abs2_cut = zeros(num_linelimits_cut)
    Yre_from_cut = zeros(num_linelimits_cut)
    Yim_from_cut = zeros(num_linelimits_cut)
    Yre_to_cut = zeros(num_linelimits_cut)
    Yim_to_cut = zeros(num_linelimits_cut)
    flowmax_cut = zeros(num_linelimits_cut)

    for i in 1:num_linelimits_cut
        # Apparent power limits (from bus)
        l = limind_cut[i]
        flowmax_cut[i] = (line_cut[l].rateA / baseMVA)^2

        Yff_abs2_cut[i] = yline_cut[l].YffR^2 + yline_cut[l].YffI^2
        Yft_abs2_cut[i] = yline_cut[l].YftR^2 + yline_cut[l].YftI^2
        Yre_from_cut[i] = yline_cut[l].YffR*yline_cut[l].YftR + yline_cut[l].YffI*yline_cut[l].YftI
        Yim_from_cut[i] = -yline_cut[l].YffR*yline_cut[l].YftI + yline_cut[l].YffI*yline_cut[l].YftR

        Ytf_abs2_cut[i] = yline_cut[l].YtfR^2 + yline_cut[l].YtfI^2
        Ytt_abs2_cut[i] = yline_cut[l].YttR^2 + yline_cut[l].YttI^2
        Yre_to_cut[i] = yline_cut[l].YtfR*yline_cut[l].YttR + yline_cut[l].YtfI*yline_cut[l].YttI
        Yim_to_cut[i] = -yline_cut[l].YtfR*yline_cut[l].YttI + yline_cut[l].YtfI*yline_cut[l].YttR
    end
    # Constraints

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])
    @variable(m, y_from[t=1:cut_time, i=1:num_linelimits] >= 0)
    @variable(m, y_to[t=1:cut_time, i=1:num_linelimits] >= 0)
    @variable(m, y_from_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)
    @variable(m, y_to_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            @constraint(m, ramping[t=1:T-1,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    #= Objective function
    =#
#    f(x, y) = max(x, y)
#    JuMP.register(m, :f, 2, f, autodiff=true)
    @NLobjective(m, Min,
         sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2
            + gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g])
            + gen[g].coeff[gen[g].n] for t=1:T,g=1:num_gens)  # FIXME
        + lambda_p * sum(  # real power balance until cut
          ((sum(yline[l].YffR for l in frombus[b])
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
          - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA)^2
          for t=1:cut_time, b=1:num_buses
        )
        + lambda_q * sum(  # reactive power balance until cut
          ((sum(-yline[l].YffI for l in frombus[b])
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
          - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA)^2
          for t=1:cut_time, b=1:num_buses
        )
        + lambda_p * sum(  # real power balance after cut
          ((sum(yline_cut[l].YffR for l in frombus_cut[b])
           + sum(yline_cut[l].YttR for l in tobus_cut[b])
           + ybus[b].YshR)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (yline_cut[l].YftR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (yline_cut[l].YtfR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA)^2
          for t=cut_time+1:T, b=1:num_buses
        )
        + lambda_q * sum(  # reactive power balance after cut
          ((sum(-yline_cut[l].YffI for l in frombus_cut[b])
           + sum(-yline_cut[l].YttI for l in tobus_cut[b])
           - ybus[b].YshI)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (-yline_cut[l].YftI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (-yline_cut[l].YtfI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA)^2
          for t=cut_time+1:T, b=1:num_buses
        )
        + lambda_l * sum(  # linelimit to before cut
             y_to[t, i] for t=1:cut_time, i=1:num_linelimits
        )
        + lambda_l * sum(  # linelimit from before cut
             y_from[t, i] for t=1:cut_time, i=1:num_linelimits
        )
        + lambda_l * sum(  # linelimit to after cut
             y_to_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
        + lambda_l * sum(  # linelimit from after cut
             y_from_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
    )
#        + T * lambda_l * sum(  # linelimit to until cut
#          f(0,
#              Vm[t,busdict[line[limind[i]].from]]^2 *
#              (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
#               + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
#               + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
#               (Yre_to[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
#                - Yim_to[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
#               ) - flowmax[i]
#          )^2
#          for t=1:cut_time, i=1:num_linelimits
#        )
#        + T * lambda_l * sum(  # linelimit from until cut
#          f(0,
#              Vm[t,busdict[line[limind[i]].to]]^2 *
#              (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
#               + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
#               + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
#               (Yre_from[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
#                -Yim_from[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
#               ) - flowmax[i]
#          )^2
#          for t=1:cut_time, i=1:num_linelimits
#        )
#        + T * lambda_l * sum(  # linelimit to after cut
#          f(0,
#              Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2 *
#              (Ytf_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
#               + Ytt_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
#               + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
#               (Yre_to_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
#                - Yim_to_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
#               ) - flowmax_cut[i]
#          )^2
#          for t=cut_time+1:T, i=1:num_linelimits_cut
#        )
#        + T * lambda_l * sum(  # linelimit from after cut
#          f(0,
#              Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2 *
#              (Yff_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
#               + Yft_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
#               + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
#               (Yre_from_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
#                -Yim_from_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
#               ) - flowmax_cut[i]
#          )^2
#          for t=cut_time+1:T, i=1:num_linelimits_cut
#        )

    # linelimit to until cut
    @NLconstraint(m, limit_to[t=1:cut_time, i=1:num_linelimits],
        y_to[t, i] >=
        Vm[t,busdict[line[limind[i]].from]]^2 *
        (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
        + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
        + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
        (Yre_to[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
        - Yim_to[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
        ) - flowmax[i]
    )
    @NLconstraint(m, limit_from[t=1:cut_time, i=1:num_linelimits],
        y_from[t, i] >=
        Vm[t,busdict[line[limind[i]].to]]^2 *
        (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
        + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
        + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
        (Yre_from[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
        -Yim_from[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
        ) - flowmax[i]
      )
    #linelimit to after cut
    @NLconstraint(m, cut_limit_to[t=cut_time+1:T, i=1:num_linelimits_cut],
        y_to_cut[t, i] >=
        Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2 *
        (Ytf_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
        + Ytt_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
        + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
        (Yre_to_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
        - Yim_to_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
        ) - flowmax_cut[i]
    )
    #linelimit from after cut
    @NLconstraint(m, cut_limit_from[t=cut_time+1:T, i=1:num_linelimits_cut],
        y_from_cut[t, i] >=
        Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2 *
        (Yff_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
        + Yft_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
        + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
        (Yre_from_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
        -Yim_from_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
        ) - flowmax_cut[i]
    )
    return m
end

function get_cut_mpcmodel2(circuit, cut_circuit, demand, cut_time; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing)
    # Model2 Hard constraint until cut_time
    # Soft constraint after cut_time
    # No Power generation cost
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    lambda_p = 3000
    lambda_q = 3000
    lambda_l = 3000

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
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)

    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])
    @variable(m, y_from_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)
    @variable(m, y_to_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            @constraint(m, ramping[t=1:T-1,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    # Line limits
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)
    Yre_from = zeros(num_linelimits)
    Yim_from = zeros(num_linelimits)
    Yre_to = zeros(num_linelimits)
    Yim_to = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2

        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre_from[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim_from[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR

        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre_to[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim_to[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end


    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)

    Yff_abs2_cut = zeros(num_linelimits_cut)
    Yft_abs2_cut = zeros(num_linelimits_cut)
    Ytf_abs2_cut = zeros(num_linelimits_cut)
    Ytt_abs2_cut = zeros(num_linelimits_cut)
    Yre_from_cut = zeros(num_linelimits_cut)
    Yim_from_cut = zeros(num_linelimits_cut)
    Yre_to_cut = zeros(num_linelimits_cut)
    Yim_to_cut = zeros(num_linelimits_cut)
    flowmax_cut = zeros(num_linelimits_cut)

    for i in 1:num_linelimits_cut
        # Apparent power limits (from bus)
        l = limind_cut[i]
        flowmax_cut[i] = (line_cut[l].rateA / baseMVA)^2

        Yff_abs2_cut[i] = yline_cut[l].YffR^2 + yline_cut[l].YffI^2
        Yft_abs2_cut[i] = yline_cut[l].YftR^2 + yline_cut[l].YftI^2
        Yre_from_cut[i] = yline_cut[l].YffR*yline_cut[l].YftR + yline_cut[l].YffI*yline_cut[l].YftI
        Yim_from_cut[i] = -yline_cut[l].YffR*yline_cut[l].YftI + yline_cut[l].YffI*yline_cut[l].YftR

        Ytf_abs2_cut[i] = yline_cut[l].YtfR^2 + yline_cut[l].YtfI^2
        Ytt_abs2_cut[i] = yline_cut[l].YttR^2 + yline_cut[l].YttI^2
        Yre_to_cut[i] = yline_cut[l].YtfR*yline_cut[l].YttR + yline_cut[l].YtfI*yline_cut[l].YttI
        Yim_to_cut[i] = -yline_cut[l].YtfR*yline_cut[l].YttI + yline_cut[l].YtfI*yline_cut[l].YttR
    end
    #= Objective function
    =#
#    f(x, y) = max(x, y)
#    JuMP.register(m, :f, 2, f, autodiff=true)  # Then you can use f(x, 0) = max(x, 0)
    @NLobjective(m, Min,
#         sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2
#            + gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g])
#            + gen[g].coeff[gen[g].n] for t=1:cut_time,g=1:num_gens)  # FIXME
        lambda_p * sum(  # real power balance after cut
          ((sum(yline_cut[l].YffR for l in frombus_cut[b])
           + sum(yline_cut[l].YttR for l in tobus_cut[b])
           + ybus[b].YshR)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (yline_cut[l].YftR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (yline_cut[l].YtfR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA)^2
          for t=cut_time+1:T, b=1:num_buses
        )
        + lambda_q * sum(  # reactive power balance after cut
          ((sum(-yline_cut[l].YffI for l in frombus_cut[b])
           + sum(-yline_cut[l].YttI for l in tobus_cut[b])
           - ybus[b].YshI)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (-yline_cut[l].YftI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (-yline_cut[l].YtfI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA)^2
          for t=cut_time+1:T, b=1:num_buses
        )
        + lambda_l * sum(  # linelimit to after cut
             y_to_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
        + lambda_l * sum(  # linelimit from after cut
             y_from_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
    )

    # Power flow constraints before cut: real part
    @NLconstraint(m, pfreal[t=1:cut_time,b=1:num_buses],
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

    # Power flow constraints before cut: imaginary part
    @NLconstraint(m, pfimag[t=1:cut_time,b=1:num_buses],
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

    # Linelimit from before cut
    @NLconstraint(m, flowmaxfrom[t=1:cut_time,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].from]]^2 *
                  (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre_from[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    - Yim_from[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <= 0)

    # Linelimit to before cut
    @NLconstraint(m, flowmaxto[t=1:cut_time,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].to]]^2 *
                  (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre_to[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    -Yim_to[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <=0)
    #linelimit from after cut
    @NLconstraint(m, cut_limit_from[t=cut_time+1:T, i=1:num_linelimits_cut],
          y_from_cut[t, i] >=
          Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2 *
          (Yff_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
           + Yft_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
           + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
           (Yre_from_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
            -Yim_from_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
           ) - flowmax_cut[i]
    )
    #linelimit to after cut
    @NLconstraint(m, cut_limit_to[t=cut_time+1:T, i=1:num_linelimits_cut],
          y_to_cut[t, i] >=
          Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2 *
          (Ytf_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
           + Ytt_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
           + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
           (Yre_to_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
            - Yim_to_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
           ) - flowmax_cut[i]
    )
    return m
end


function get_cut_mpcmodel3(circuit, cut_circuit, demand, cut_time; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing)
    # Model3 Soft constraint before cut and after cut
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    lambda_p = 3000
    lambda_q = 3000
    lambda_l = 3000

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
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)

    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])
    @variable(m, y_from_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)
    @variable(m, y_to_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            @constraint(m, ramping[t=1:T-1,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    # Line limits
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)
    Yre_from = zeros(num_linelimits)
    Yim_from = zeros(num_linelimits)
    Yre_to = zeros(num_linelimits)
    Yim_to = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2

        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre_from[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim_from[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR

        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre_to[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim_to[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end


    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)

    Yff_abs2_cut = zeros(num_linelimits_cut)
    Yft_abs2_cut = zeros(num_linelimits_cut)
    Ytf_abs2_cut = zeros(num_linelimits_cut)
    Ytt_abs2_cut = zeros(num_linelimits_cut)
    Yre_from_cut = zeros(num_linelimits_cut)
    Yim_from_cut = zeros(num_linelimits_cut)
    Yre_to_cut = zeros(num_linelimits_cut)
    Yim_to_cut = zeros(num_linelimits_cut)
    flowmax_cut = zeros(num_linelimits_cut)

    for i in 1:num_linelimits_cut
        # Apparent power limits (from bus)
        l = limind_cut[i]
        flowmax_cut[i] = (line_cut[l].rateA / baseMVA)^2

        Yff_abs2_cut[i] = yline_cut[l].YffR^2 + yline_cut[l].YffI^2
        Yft_abs2_cut[i] = yline_cut[l].YftR^2 + yline_cut[l].YftI^2
        Yre_from_cut[i] = yline_cut[l].YffR*yline_cut[l].YftR + yline_cut[l].YffI*yline_cut[l].YftI
        Yim_from_cut[i] = -yline_cut[l].YffR*yline_cut[l].YftI + yline_cut[l].YffI*yline_cut[l].YftR

        Ytf_abs2_cut[i] = yline_cut[l].YtfR^2 + yline_cut[l].YtfI^2
        Ytt_abs2_cut[i] = yline_cut[l].YttR^2 + yline_cut[l].YttI^2
        Yre_to_cut[i] = yline_cut[l].YtfR*yline_cut[l].YttR + yline_cut[l].YtfI*yline_cut[l].YttI
        Yim_to_cut[i] = -yline_cut[l].YtfR*yline_cut[l].YttI + yline_cut[l].YtfI*yline_cut[l].YttR
    end
    #= Objective function
    =#
    f(x, y) = max(x, y)
    JuMP.register(m, :f, 2, f, autodiff=true)
    @NLobjective(m, Min,
         sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2
            + gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g])
            + gen[g].coeff[gen[g].n] for t=1:cut_time,g=1:num_gens)  # FIXME
        + lambda_l * sum(  # linelimit to after cut
             y_to_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
        + lambda_l * sum(  # linelimit from after cut
             y_from_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
    )

    # before cut Power flow constraints: real part
    @NLconstraint(m, pfreal[t=1:cut_time,b=1:num_buses],
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
    @NLconstraint(m, pfimag[t=1:cut_time,b=1:num_buses],
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

    # after cut power flow constraints
    @NLconstraint(m, pfreal_cut[t=cut_time+1:T,b=1:num_buses],
        ((sum(yline_cut[l].YffR for l in frombus_cut[b])
           + sum(yline_cut[l].YttR for l in tobus_cut[b])
           + ybus[b].YshR)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (yline_cut[l].YftR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (yline_cut[l].YtfR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA
        ) == 0
    )

    @NLconstraint(m, pfimag_cut[t=cut_time+1:T,b=1:num_buses],
        ((sum(-yline_cut[l].YffI for l in frombus_cut[b])
           + sum(-yline_cut[l].YttI for l in tobus_cut[b])
           - ybus[b].YshI)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (-yline_cut[l].YftI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (-yline_cut[l].YtfI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA
        ) == 0
    )

    # linelimit before cut
    @NLconstraint(m, flowmaxfrom[t=1:cut_time,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].from]]^2 *
                  (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre_from[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    - Yim_from[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <= 0)

    @NLconstraint(m, flowmaxto[t=1:cut_time,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].to]]^2 *
                  (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre_to[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    -Yim_to[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <=0)

    #TO after cut
    @NLconstraint(m, cut_limit_to[t=cut_time+1:T, i=1:num_linelimits_cut],
          y_to_cut[t, i] >=
          Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2 *
          (Ytf_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
           + Ytt_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
           + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
           (Yre_to_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
            - Yim_to_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
           ) - flowmax_cut[i]
    )
    #From after cut
    @NLconstraint(m, cut_limit_from[t=cut_time+1:T, i=1:num_linelimits_cut],
          y_from_cut[t, i] >=
          Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2 *
          (Yff_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
           + Yft_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
           + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
           (Yre_from_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
            -Yim_from_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
           ) - flowmax_cut[i]
    )
    return m
end


function get_cut_mpcmodel4(circuit, cut_circuit, demand, cut_time; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing)
    # Model4 Soft constraints for powerflow; no objective; no line constraints
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    # TODO: Need to receive cut circuit and cut time and divide it with before cut and after cut
    lambda_p = 3000
    lambda_q = 3000
#    lambda_l = 3000

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
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)

    # Line limits
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)
    Yre_from = zeros(num_linelimits)
    Yim_from = zeros(num_linelimits)
    Yre_to = zeros(num_linelimits)
    Yim_to = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2

        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre_from[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim_from[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR

        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre_to[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim_to[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end


    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)

    Yff_abs2_cut = zeros(num_linelimits_cut)
    Yft_abs2_cut = zeros(num_linelimits_cut)
    Ytf_abs2_cut = zeros(num_linelimits_cut)
    Ytt_abs2_cut = zeros(num_linelimits_cut)
    Yre_from_cut = zeros(num_linelimits_cut)
    Yim_from_cut = zeros(num_linelimits_cut)
    Yre_to_cut = zeros(num_linelimits_cut)
    Yim_to_cut = zeros(num_linelimits_cut)
    flowmax_cut = zeros(num_linelimits_cut)

    for i in 1:num_linelimits_cut
        # Apparent power limits (from bus)
        l = limind_cut[i]
        flowmax_cut[i] = (line_cut[l].rateA / baseMVA)^2

        Yff_abs2_cut[i] = yline_cut[l].YffR^2 + yline_cut[l].YffI^2
        Yft_abs2_cut[i] = yline_cut[l].YftR^2 + yline_cut[l].YftI^2
        Yre_from_cut[i] = yline_cut[l].YffR*yline_cut[l].YftR + yline_cut[l].YffI*yline_cut[l].YftI
        Yim_from_cut[i] = -yline_cut[l].YffR*yline_cut[l].YftI + yline_cut[l].YffI*yline_cut[l].YftR

        Ytf_abs2_cut[i] = yline_cut[l].YtfR^2 + yline_cut[l].YtfI^2
        Ytt_abs2_cut[i] = yline_cut[l].YttR^2 + yline_cut[l].YttI^2
        Yre_to_cut[i] = yline_cut[l].YtfR*yline_cut[l].YttR + yline_cut[l].YtfI*yline_cut[l].YttI
        Yim_to_cut[i] = -yline_cut[l].YtfR*yline_cut[l].YttI + yline_cut[l].YtfI*yline_cut[l].YttR
    end
    # Constraints

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])
#    @variable(m, y_from[t=1:cut_time, i=1:num_linelimits] >= 0)
#    @variable(m, y_to[t=1:cut_time, i=1:num_linelimits] >= 0)
#    @variable(m, y_from_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)
#    @variable(m, y_to_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            @constraint(m, ramping[t=1:T-1,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    #= Objective function
    =#
    @NLobjective(m, Min,
        lambda_p * sum(  # real power balance until cut
          ((sum(yline[l].YffR for l in frombus[b])
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
          - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA)^2
          for t=1:cut_time, b=1:num_buses
        )
        + lambda_q * sum(  # reactive power balance until cut
          ((sum(-yline[l].YffI for l in frombus[b])
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
          - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA)^2
          for t=1:cut_time, b=1:num_buses
        )
        + lambda_p * sum(  # real power balance after cut
          ((sum(yline_cut[l].YffR for l in frombus_cut[b])
           + sum(yline_cut[l].YttR for l in tobus_cut[b])
           + ybus[b].YshR)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (yline_cut[l].YftR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (yline_cut[l].YtfR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA)^2
          for t=cut_time+1:T, b=1:num_buses
        )
        + lambda_q * sum(  # reactive power balance after cut
          ((sum(-yline_cut[l].YffI for l in frombus_cut[b])
           + sum(-yline_cut[l].YttI for l in tobus_cut[b])
           - ybus[b].YshI)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (-yline_cut[l].YftI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (-yline_cut[l].YtfI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA)^2
          for t=cut_time+1:T, b=1:num_buses
        )
    )
    return m
end


function get_cut_mpcmodel5(circuit, cut_circuit, demand, cut_time; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing)
    # Model5 Soft constraints all time, no objective function
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    lambda_p = 3000
    lambda_q = 3000
    lambda_l = 3000

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
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)

    # Line limits
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)
    Yre_from = zeros(num_linelimits)
    Yim_from = zeros(num_linelimits)
    Yre_to = zeros(num_linelimits)
    Yim_to = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2

        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre_from[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim_from[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR

        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre_to[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim_to[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end


    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)

    Yff_abs2_cut = zeros(num_linelimits_cut)
    Yft_abs2_cut = zeros(num_linelimits_cut)
    Ytf_abs2_cut = zeros(num_linelimits_cut)
    Ytt_abs2_cut = zeros(num_linelimits_cut)
    Yre_from_cut = zeros(num_linelimits_cut)
    Yim_from_cut = zeros(num_linelimits_cut)
    Yre_to_cut = zeros(num_linelimits_cut)
    Yim_to_cut = zeros(num_linelimits_cut)
    flowmax_cut = zeros(num_linelimits_cut)

    for i in 1:num_linelimits_cut
        # Apparent power limits (from bus)
        l = limind_cut[i]
        flowmax_cut[i] = (line_cut[l].rateA / baseMVA)^2

        Yff_abs2_cut[i] = yline_cut[l].YffR^2 + yline_cut[l].YffI^2
        Yft_abs2_cut[i] = yline_cut[l].YftR^2 + yline_cut[l].YftI^2
        Yre_from_cut[i] = yline_cut[l].YffR*yline_cut[l].YftR + yline_cut[l].YffI*yline_cut[l].YftI
        Yim_from_cut[i] = -yline_cut[l].YffR*yline_cut[l].YftI + yline_cut[l].YffI*yline_cut[l].YftR

        Ytf_abs2_cut[i] = yline_cut[l].YtfR^2 + yline_cut[l].YtfI^2
        Ytt_abs2_cut[i] = yline_cut[l].YttR^2 + yline_cut[l].YttI^2
        Yre_to_cut[i] = yline_cut[l].YtfR*yline_cut[l].YttR + yline_cut[l].YtfI*yline_cut[l].YttI
        Yim_to_cut[i] = -yline_cut[l].YtfR*yline_cut[l].YttI + yline_cut[l].YtfI*yline_cut[l].YttR
    end
    # Constraints

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])
    @variable(m, y_from[t=1:cut_time, i=1:num_linelimits] >= 0)
    @variable(m, y_to[t=1:cut_time, i=1:num_linelimits] >= 0)
    @variable(m, y_from_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)
    @variable(m, y_to_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            @constraint(m, ramping[t=1:T-1,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    #= Objective function
    =#
#    f(x, y) = max(x, y)
#    JuMP.register(m, :f, 2, f, autodiff=true)
    @NLobjective(m, Min,
        lambda_p * sum(  # real power balance until cut
          ((sum(yline[l].YffR for l in frombus[b])
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
          - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA)^2
          for t=1:cut_time, b=1:num_buses
        )
        + lambda_q * sum(  # reactive power balance until cut
          ((sum(-yline[l].YffI for l in frombus[b])
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
          - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA)^2
          for t=1:cut_time, b=1:num_buses
        )
        + lambda_p * sum(  # real power balance after cut
          ((sum(yline_cut[l].YffR for l in frombus_cut[b])
           + sum(yline_cut[l].YttR for l in tobus_cut[b])
           + ybus[b].YshR)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (yline_cut[l].YftR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (yline_cut[l].YtfR*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfI*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA)^2
          for t=cut_time+1:T, b=1:num_buses
        )
        + lambda_q * sum(  # reactive power balance after cut
          ((sum(-yline_cut[l].YffI for l in frombus_cut[b])
           + sum(-yline_cut[l].YttI for l in tobus_cut[b])
           - ybus[b].YshI)*Vm[t,b]^2
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].to]]*
                (-yline_cut[l].YftI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]])
                 + yline_cut[l].YftR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].to]]))
                for l in frombus_cut[b])
          + sum(Vm[t,b]*Vm[t,busdict_cut[line_cut[l].from]]*
                (-yline_cut[l].YtfI*cos(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]])
                 + yline_cut[l].YtfR*sin(Va[t,b]-Va[t,busdict_cut[line_cut[l].from]]))
                for l in tobus_cut[b])
          - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA)^2
          for t=cut_time+1:T, b=1:num_buses
        )
        + lambda_l * sum(  # linelimit to before cut
             y_to[t, i] for t=1:cut_time, i=1:num_linelimits
        )
        + lambda_l * sum(  # linelimit from before cut
             y_from[t, i] for t=1:cut_time, i=1:num_linelimits
        )
        + lambda_l * sum(  # linelimit to after cut
             y_to_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
        + lambda_l * sum(  # linelimit from after cut
             y_from_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
    )
    # linelimit to until cut
    @NLconstraint(m, limit_to[t=1:cut_time, i=1:num_linelimits],
        y_to[t, i] >=
        Vm[t,busdict[line[limind[i]].from]]^2 *
        (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
        + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
        + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
        (Yre_to[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
        - Yim_to[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
        ) - flowmax[i]
    )
    # linelimit from until cut
    @NLconstraint(m, limit_from[t=1:cut_time, i=1:num_linelimits],
        y_from[t, i] >=
        Vm[t,busdict[line[limind[i]].to]]^2 *
        (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
        + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
        + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
        (Yre_from[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
        -Yim_from[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
        ) - flowmax[i]
      )
    #linelimit to after cut
    @NLconstraint(m, cut_limit_to[t=cut_time+1:T, i=1:num_linelimits_cut],
        y_to_cut[t, i] >=
        Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2 *
        (Ytf_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
        + Ytt_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
        + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
        (Yre_to_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
        - Yim_to_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
        ) - flowmax_cut[i]
    )
    #linelimit from after cut
    @NLconstraint(m, cut_limit_from[t=cut_time+1:T, i=1:num_linelimits_cut],
        y_from_cut[t, i] >=
        Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2 *
        (Yff_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
        + Yft_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
        + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
        (Yre_from_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
        -Yim_from_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
        ) - flowmax_cut[i]
    )
    return m
end


function get_cut_mpcmodel6(circuit, cut_circuit, demand, cut_time; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing)
    # Model6 Soft constraints all time, no objective function and no powerflow constraints
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    lambda_p = 3000
    lambda_q = 3000
    lambda_l = 3000

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
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)

    # Line limits
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)
    Yre_from = zeros(num_linelimits)
    Yim_from = zeros(num_linelimits)
    Yre_to = zeros(num_linelimits)
    Yim_to = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2

        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre_from[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim_from[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR

        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre_to[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim_to[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end


    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)

    Yff_abs2_cut = zeros(num_linelimits_cut)
    Yft_abs2_cut = zeros(num_linelimits_cut)
    Ytf_abs2_cut = zeros(num_linelimits_cut)
    Ytt_abs2_cut = zeros(num_linelimits_cut)
    Yre_from_cut = zeros(num_linelimits_cut)
    Yim_from_cut = zeros(num_linelimits_cut)
    Yre_to_cut = zeros(num_linelimits_cut)
    Yim_to_cut = zeros(num_linelimits_cut)
    flowmax_cut = zeros(num_linelimits_cut)

    for i in 1:num_linelimits_cut
        # Apparent power limits (from bus)
        l = limind_cut[i]
        flowmax_cut[i] = (line_cut[l].rateA / baseMVA)^2

        Yff_abs2_cut[i] = yline_cut[l].YffR^2 + yline_cut[l].YffI^2
        Yft_abs2_cut[i] = yline_cut[l].YftR^2 + yline_cut[l].YftI^2
        Yre_from_cut[i] = yline_cut[l].YffR*yline_cut[l].YftR + yline_cut[l].YffI*yline_cut[l].YftI
        Yim_from_cut[i] = -yline_cut[l].YffR*yline_cut[l].YftI + yline_cut[l].YffI*yline_cut[l].YftR

        Ytf_abs2_cut[i] = yline_cut[l].YtfR^2 + yline_cut[l].YtfI^2
        Ytt_abs2_cut[i] = yline_cut[l].YttR^2 + yline_cut[l].YttI^2
        Yre_to_cut[i] = yline_cut[l].YtfR*yline_cut[l].YttR + yline_cut[l].YtfI*yline_cut[l].YttI
        Yim_to_cut[i] = -yline_cut[l].YtfR*yline_cut[l].YttI + yline_cut[l].YtfI*yline_cut[l].YttR
    end
    # Constraints

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])
    @variable(m, y_from[t=1:cut_time, i=1:num_linelimits] >= 0)
    @variable(m, y_to[t=1:cut_time, i=1:num_linelimits] >= 0)
    @variable(m, y_from_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)
    @variable(m, y_to_cut[t=cut_time+1:T, i=1:num_linelimits_cut] >= 0)

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            @constraint(m, ramping[t=1:T-1,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    #= Objective function
    =#
#    f(x, y) = max(x, y)
#    JuMP.register(m, :f, 2, f, autodiff=true)
    @NLobjective(m, Min,
        lambda_l * sum(  # linelimit to before cut
             y_to[t, i] for t=1:cut_time, i=1:num_linelimits
        )
        + lambda_l * sum(  # linelimit from before cut
             y_from[t, i] for t=1:cut_time, i=1:num_linelimits
        )
        + lambda_l * sum(  # linelimit to after cut
             y_to_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
        + lambda_l * sum(  # linelimit from after cut
             y_from_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
        )
    )
    # linelimit to until cut
    @NLconstraint(m, limit_to[t=1:cut_time, i=1:num_linelimits],
        y_to[t, i] >=
        Vm[t,busdict[line[limind[i]].from]]^2 *
        (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
        + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
        + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
        (Yre_to[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
        - Yim_to[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
        ) - flowmax[i]
    )
    # linelimit from until cut
    @NLconstraint(m, limit_from[t=1:cut_time, i=1:num_linelimits],
        y_from[t, i] >=
        Vm[t,busdict[line[limind[i]].to]]^2 *
        (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
        + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
        + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
        (Yre_from[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
        -Yim_from[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
        ) - flowmax[i]
      )
    #linelimit to after cut
    @NLconstraint(m, cut_limit_to[t=cut_time+1:T, i=1:num_linelimits_cut],
        y_to_cut[t, i] >=
        Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2 *
        (Ytf_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
        + Ytt_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
        + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
        (Yre_to_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
        - Yim_to_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
        ) - flowmax_cut[i]
    )
    #linelimit from after cut
    @NLconstraint(m, cut_limit_from[t=cut_time+1:T, i=1:num_linelimits_cut],
        y_from_cut[t, i] >=
        Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2 *
        (Yff_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
        + Yft_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
        + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
        (Yre_from_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
        -Yim_from_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
        ) - flowmax_cut[i]
    )
    return m
end


function get_cut_mpcmodel7(circuit, cut_circuit, demand, cut_time; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing)
    # Model7 Only hard line constraints (Try without ramping)
    #=
    Args:
        circuit - Circuit instance
        demand - Load instance
    =#

    lambda_p = 3000
    lambda_q = 3000
    lambda_l = 3000

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
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)

    # Line limits
    rateA = getfield.(line, :rateA)  # equiv to getattr for struct
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)
    Yre_from = zeros(num_linelimits)
    Yim_from = zeros(num_linelimits)
    Yre_to = zeros(num_linelimits)
    Yim_to = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2

        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre_from[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim_from[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR

        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre_to[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim_to[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end


    # Shortcuts for cut circuit
    line_cut = cut_circuit.line
    yline_cut = cut_circuit.yline
    busdict_cut = cut_circuit.busdict
    frombus_cut = cut_circuit.frombus
    tobus_cut = cut_circuit.tobus

    # Line limits
    rateA_cut = getfield.(line_cut, :rateA)  # equiv to getattr for struct
    limind_cut = findall((rateA_cut .!= 0) .& (rateA_cut .< 1.0e10))
    num_linelimits_cut = length(limind_cut)
    println("########################################################")
    println("Cut num linelimits: $(num_linelimits_cut)")
    println("Num linelimits: $(num_linelimits)")
    println("########################################################")

    Yff_abs2_cut = zeros(num_linelimits_cut)
    Yft_abs2_cut = zeros(num_linelimits_cut)
    Ytf_abs2_cut = zeros(num_linelimits_cut)
    Ytt_abs2_cut = zeros(num_linelimits_cut)
    Yre_from_cut = zeros(num_linelimits_cut)
    Yim_from_cut = zeros(num_linelimits_cut)
    Yre_to_cut = zeros(num_linelimits_cut)
    Yim_to_cut = zeros(num_linelimits_cut)
    flowmax_cut = zeros(num_linelimits_cut)

    for i in 1:num_linelimits_cut
        # Apparent power limits (from bus)
        l = limind_cut[i]
        flowmax_cut[i] = (line_cut[l].rateA / baseMVA)^2

        Yff_abs2_cut[i] = yline_cut[l].YffR^2 + yline_cut[l].YffI^2
        Yft_abs2_cut[i] = yline_cut[l].YftR^2 + yline_cut[l].YftI^2
        Yre_from_cut[i] = yline_cut[l].YffR*yline_cut[l].YftR + yline_cut[l].YffI*yline_cut[l].YftI
        Yim_from_cut[i] = -yline_cut[l].YffR*yline_cut[l].YftI + yline_cut[l].YffI*yline_cut[l].YftR

        Ytf_abs2_cut[i] = yline_cut[l].YtfR^2 + yline_cut[l].YtfI^2
        Ytt_abs2_cut[i] = yline_cut[l].YttR^2 + yline_cut[l].YttI^2
        Yre_to_cut[i] = yline_cut[l].YtfR*yline_cut[l].YttR + yline_cut[l].YtfI*yline_cut[l].YttI
        Yim_to_cut[i] = -yline_cut[l].YtfR*yline_cut[l].YttI + yline_cut[l].YtfI*yline_cut[l].YttR
    end
    # Constraints

    #=
    Pg[t, g]: real power at (time, generator) p_{t,g}
    Qg[t, g]: reactive power at (time, generator) q_{t,g}
    Vm[t, b]: voltage magnitude at (time, bus) v_{t,i}
    Va[t, b]: voltage angle at (time, bus) \theta_{t, i}
    =#
    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])

    # Voltage angle has fixed value at reference bus (lower = upper)
    # But no lower or upper bound on other Va[t, b]. Why? because it is sin, cos so doesn't matter.
    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
            @constraint(m, ramping[t=1:T-1,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
        else
            @constraint(m, ramping[t=1:T,g=1:num_gens],
                    -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    #= Objective function
    =#
#    @NLobjective(m, Min,
#        lambda_l * sum(  # linelimit to before cut
#             y_to[t, i] for t=1:cut_time, i=1:num_linelimits
#        )
#        + lambda_l * sum(  # linelimit from before cut
#             y_from[t, i] for t=1:cut_time, i=1:num_linelimits
#        )
#        + lambda_l * sum(  # linelimit to after cut
#             y_to_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
#        )
#        + lambda_l * sum(  # linelimit from after cut
#             y_from_cut[t, i] for t=cut_time+1:T, i=1:num_linelimits_cut
#        )
#    )
    # Hard Linelimit from before cut
    @NLconstraint(m, flowmaxfrom[t=1:cut_time,i=1:num_linelimits],
        Vm[t,busdict[line[limind[i]].from]]^2 *
        (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
        + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
        + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
        (Yre_from[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
        - Yim_from[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
        ) - flowmax[i] <= 0
    )

    # Hard Linelimit to before cut
    @NLconstraint(m, flowmaxto[t=1:cut_time,i=1:num_linelimits],
        Vm[t,busdict[line[limind[i]].to]]^2 *
        (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
        + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
        + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
        (Yre_to[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
        -Yim_to[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
        ) - flowmax[i] <= 0
    )
    # Hard linelimit from after cut
    @NLconstraint(m, cut_limit_from[t=cut_time+1:T, i=1:num_linelimits_cut],
        Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2 *
        (Yff_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
        + Yft_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
        + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
        (Yre_from_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
        -Yim_from_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
        ) - flowmax_cut[i] <= 0
    )
    # Hard linelimit to after cut
    @NLconstraint(m, cut_limit_to[t=cut_time+1:T, i=1:num_linelimits_cut],
        Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2 *
        (Ytf_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
        + Ytt_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
        + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
        (Yre_to_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
        - Yim_to_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
        ) - flowmax_cut[i] <= 0
    )
    return m
end


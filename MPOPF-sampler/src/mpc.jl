using FileIO
using DelimitedFiles
using JuMP
using MathProgBase
using Ipopt
using Printf
using Random
using JLD2
using Dates
using DataStructures

include("mpc_data.jl")
include("mpc_model.jl")
include("mpc_qp.jl")
include("mpc_cut_model.jl")


function runif(n, l, u)
    #=
    Args
        n (Int) - # of samples
        l (Float, Array{Float, 1}) - lower bound
        u (Float, Array{Float, 1}) - upper bound

    Returns
        u ~ U[l, u]
    =#
    return (u - l) .* rand(n) .+ l
end

function init_x(m, circuit, demand; uniform=false)
    num_buses = length(circuit.bus)
    num_gens = length(circuit.gen)
    T = size(demand.pd,2)

    Vm = zeros(T, num_buses)

    for b=1:num_buses
        Vm[:,b] .= (
            uniform ? runif(1, circuit.bus[b].Vmin, circuit.bus[b].Vmax) : 0.5*(circuit.bus[b].Vmax + circuit.bus[b].Vmin)
           )
    end
    Va = circuit.bus[circuit.busref].Va * ones(T, num_buses)

    setvalue(m[:Vm], Vm)
    setvalue(m[:Va], Va)

    Pg = zeros(T, num_gens)
    Qg = zeros(T, num_gens)

    for g=1:num_gens
        Pg[:,g] .= (
            uniform ? runif(1, circuit.gen[g].Pmin, circuit.gen[g].Pmax) : 0.5*(circuit.gen[g].Pmax + circuit.gen[g].Pmin)
           )
        Qg[:,g] .= (
            uniform ? runif(1, circuit.gen[g].Qmin, circuit.gen[g].Qmax) : 0.5*(circuit.gen[g].Qmax + circuit.gen[g].Qmin)
           )
    end

    setvalue(m[:Pg], Pg)
    setvalue(m[:Qg], Qg)
end

function init_x(m_to, m_from)
    setvalue(m_to[:Vm], getvalue(m_from[:Vm]))
    setvalue(m_to[:Va], getvalue(m_from[:Va]))

    setvalue(m_to[:Pg], getvalue(m_from[:Pg]))
    setvalue(m_to[:Qg], getvalue(m_from[:Qg]))
end

function init_x(m_to, m_from, circuit, cut_circuit, demand, cut_time, T)
    init_x(m_to, circuit, demand)
    # Constraints
    setvalue(m_to[:Vm][1:cut_time, :], getvalue(m_from[:Vm][1:cut_time, :]))
    setvalue(m_to[:Va][1:cut_time, :], getvalue(m_from[:Va][1:cut_time, :]))

    setvalue(m_to[:Pg][1:cut_time, :], getvalue(m_from[:Pg][1:cut_time, :]))
    setvalue(m_to[:Qg][1:cut_time, :], getvalue(m_from[:Qg][1:cut_time, :]))

#    setvalue(m_to[:Vm], getvalue(m_from[:Vm]))
#    setvalue(m_to[:Va], getvalue(m_from[:Va]))
#
#    setvalue(m_to[:Pg], getvalue(m_from[:Pg]))
#    setvalue(m_to[:Qg], getvalue(m_from[:Qg]))

    # Set slack variables
    baseMVA = cut_circuit.baseMVA
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

    Vm = getvalue(m_to[:Vm])
    Va = getvalue(m_to[:Va])
    y_from_cut_init = zeros(T-cut_time, num_linelimits_cut)
    y_to_cut_init = zeros(T-cut_time, num_linelimits_cut)
    for t in cut_time+1:T
        for i in 1:num_linelimits_cut
            y_from_cut_init[t-cut_time, i] = Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2 *
                (Yff_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
                + Yft_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
                + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
                (Yre_from_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
                -Yim_from_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
                ) - flowmax_cut[i]
            y_to_cut_init[t-cut_time, i] = Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2 *
                (Ytf_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]^2
                + Ytt_abs2_cut[i]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]^2
                + 2*Vm[t,busdict_cut[line_cut[limind_cut[i]].from]]*Vm[t,busdict_cut[line_cut[limind_cut[i]].to]]*
                (Yre_to_cut[i]*cos(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]])
                - Yim_to_cut[i]*sin(Va[t,busdict_cut[line_cut[limind_cut[i]].from]] - Va[t,busdict_cut[line_cut[limind_cut[i]].to]]))
               ) - flowmax_cut[i]
        end
    end

#c = 1e-1
#    setvalue(m_to[:y_from_cut][cut_time+1:T, 1:num_linelimits_cut], max.(0, y_from_cut_init) .+ c)
#    setvalue(m_to[:y_to_cut][cut_time+1:T, 1:num_linelimits_cut], max.(0, y_to_cut_init) .+ c)
    setvalue(m_to[:y_from_cut][cut_time+1:T, 1:num_linelimits_cut], y_from_cut_init)
    setvalue(m_to[:y_to_cut][cut_time+1:T, 1:num_linelimits_cut], y_to_cut_init)

end

function fix_vars_until_cut(m, circuit, cut_time)
    num_gens = length(circuit.gen)
    num_buses = length(circuit.bus)

    for i in 1:cut_time
        for j in 1:num_gens
            JuMP.fix(m[:Pg][i, j], getvalue(m[:Pg])[i, j])
            JuMP.fix(m[:Qg][i, j], getvalue(m[:Qg])[i, j])
        end
        for j in 1:num_buses
            JuMP.fix(m[:Vm][i, j], getvalue(m[:Vm])[i, j])
            JuMP.fix(m[:Va][i, j], getvalue(m[:Va])[i, j])
        end
    end
end

function get_powerflow_violation(m, circuit, demand; periods=collect(1:2))
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

    Vm = getvalue(m[:Vm])
    Va = getvalue(m[:Va])
    Pg = getvalue(m[:Pg])
    Qg = getvalue(m[:Qg])
    pfreal = zeros(length(periods), num_buses)
    pfimag = zeros(length(periods), num_buses)
    for t in periods
        for b in 1:num_buses
            frombus_b = length(frombus[b]) > 0 ? frombus[b] : 0
            tobus_b = length(tobus[b]) > 0 ? tobus[b] : 0
            bus2gen_b = length(bus2gen[b]) > 0 ? bus2gen[b] : 0

            pfreal[t-periods[1]+1, b] = (
                (sum(l > 0 ? yline[l].YffR : 0 for l in frombus_b)
                   + sum(l > 0 ? yline[l].YttR : 0 for l in tobus_b)
                   + ybus[b].YshR)*Vm[t,b]^2
                  + sum(l > 0 ? Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]])) : 0
                        for l in frombus_b)
                  + sum(l > 0 ? Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfI*sin(Va[t,b]-Va[t,busdict[line[l].from]])) : 0
                        for l in tobus_b)
                  - (sum(g > 0 ? baseMVA*Pg[t,g] : 0 for g in bus2gen_b) - Pd[b,t]) / baseMVA
            )
            pfimag[t-periods[1]+1, b] = (
                  (sum(l > 0 ? -yline[l].YffI : 0 for l in frombus_b)
                   + sum(l > 0 ? -yline[l].YttI : 0 for l in tobus_b)
                   - ybus[b].YshI)*Vm[t,b]^2
                  + sum(l > 0 ? Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (-yline[l].YftI*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftR*sin(Va[t,b]-Va[t,busdict[line[l].to]])) : 0
                        for l in frombus_b)
                  + sum(l > 0 ? Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (-yline[l].YtfI*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfR*sin(Va[t,b]-Va[t,busdict[line[l].from]])) : 0
                        for l in tobus_b)
                  - (sum(g > 0 ? baseMVA*Qg[t,g] : 0 for g in bus2gen_b) - Qd[b,t]) / baseMVA
            )
        end
    end
    return pfreal, pfimag
end

function get_profname(case, T, load_scale, ramp_scale, warm)
    profname = "ipopt_"*splitext(basename(case))[1]
    profname = profname*"_t"*string(T)
    profname = profname*"_ls"*string(load_scale)*"_rs"*string(ramp_scale)

    if warm == "cold"
        profname = profname*"_cold_phase0"
    elseif warm == "shift_copy"
        profname = profname*"_warm_phase0"
    else
        profname = profname*"_warm_phase1"
    end

    timenow = Dates.format(Dates.now(), "yyyymmddHHMMSSsss")
    profname = profname * "_$(timenow)"

    return profname
end

function copyvars(inner_to, to_start, inner_from, from_start, n)
    copyto!(inner_to.x, to_start, inner_from.x, from_start, n)
    copyto!(inner_to.mult_x_L, to_start, inner_from.mult_x_L, from_start, n)
    copyto!(inner_to.mult_x_U, to_start, inner_from.mult_x_U, from_start, n)
end

function copyconstrs(inner_to, to_start, inner_from, from_start, n)
    copyto!(inner_to.g, to_start, inner_from.g, from_start, n)
    copyto!(inner_to.mult_g, to_start, inner_from.mult_g, from_start, n)
end

function copy_modelvars(m_to, m_from, circuit, demand)
    num_gens = length(circuit.gen)
    num_buses = length(circuit.bus)
    T = size(demand.pd, 2)

    inner_to = internalmodel(m_to).inner
    inner_from = internalmodel(m_from).inner

    start_to = linearindex(m_to[:Pg][1])
    start_from = linearindex(m_from[:Pg][1])
    copyvars(inner_to, start_to, inner_from, start_from, num_gens*T)

    start_to = linearindex(m_to[:Qg][1])
    start_from = linearindex(m_from[:Qg][1])
    copyvars(inner_to, start_to, inner_from, start_from, num_gens*T)

    start_to = linearindex(m_to[:Vm][1])
    start_from = linearindex(m_from[:Vm][1])
    copyvars(inner_to, start_to, inner_from, start_from, num_buses*T)

    start_to = linearindex(m_to[:Va][1])
    start_from = linearindex(m_from[:Va][1])
    copyvars(inner_to, start_to, inner_from, start_from, num_buses*T)
end

function get_cutgen(m, circuit, T)
    Pg = getvalue(m[:Pg])

    if length(circuit.bus) == 300 || length(circuit.bus) == 73
	active = findall(x -> x > 0 && x <= 1, Pg[1,:])
    else
	active = findall(x -> x > 0, Pg[1,:])
    end

    max_g = argmax(Pg[1,active])
    min_g = argmin(Pg[1,active])
    @printf("Generataor (max) %d: %.2f MW\n", active[max_g], Pg[1,active[max_g]]*circuit.baseMVA)
    @printf("Generataor (min) %d: %.2f MW\n", active[min_g], Pg[1,active[min_g]]*circuit.baseMVA)

    return active[max_g]
end

function get_cutline(m, circuit, T)
    num_buses = length(circuit.bus)

    line = circuit.line
    ybus = circuit.ybus
    yline = circuit.yline
    busdict = circuit.busdict
    frombus = circuit.frombus
    tobus = circuit.tobus

    fromflow = zeros(T, length(line))
    toflow = zeros(T, length(line))

    Vm = getvalue(m[:Vm])
    Va = getvalue(m[:Va])

    for t=1:T,b=1:num_buses
        for l in frombus[b]
            fromflow[t,l] = (yline[l].YffR + ybus[b].YshR)*Vm[t,b]^2 +
                Vm[t,b]*Vm[t,busdict[line[l].to]]*
            (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]]) +
             yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
        end

        for l in tobus[b]
            toflow[t,l] = (yline[l].YttR + ybus[b].YshR)*Vm[t,b]^2 +
                Vm[t,b]*Vm[t,busdict[line[l].from]]*
            (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]]) +
             yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
        end
    end

    abs_fromflow = abs.(fromflow)
    idx1 = findall(abs_fromflow .> 0)
    r1 = argmin(abs_fromflow[idx1])
    @printf("Line (%5d,%5d) has absolute minimum from-flow %10.8f\n",
            line[idx1[r1][2]].from, line[idx1[r1][2]].to,
            abs_fromflow[idx1[r1]])

    abs_toflow = abs.(toflow)
    idx2 = findall(abs_toflow .> 0)
    r2 = argmin(abs_toflow[idx2])
    @printf("Line (%5d,%5d) has absolute minimum to-flow %10.8f\n",
            line[idx2[r2][2]].from, line[idx2[r2][2]].to,
            abs_toflow[idx2[r2]])

    if abs_fromflow[idx1[r1]] < abs_toflow[idx2[r2]]
        return idx1[r1][2]
    else
        return idx2[r2][2]
    end
end

function set_warm(m_cur, m_prev, circuit, demand; warm_type=:shift_copy,
                  profname="noname", piecewise=false)
    num_buses = length(circuit.bus)
    num_gens = length(circuit.gen)
    num_linconstr = MathProgBase.numlinconstr(m_cur)

    rateA = getfield.(circuit.line, :rateA)
    num_linelimits = length(findall((rateA .!= 0) .& (rateA .< 1.0e10)))

    T = size(demand.pd, 2)
    @assert T >= 2

    inner_cur = internalmodel(m_cur).inner
    inner_prev = internalmodel(m_prev).inner

    size_pl = 0
    if piecewise
        # They all have the same number of pieces so we use the first one.
        size_pl = num_gens*(circuit.gen[1].n - 1)
    end

    # -----------------------------------------------------------------
    # Copy variables and their multipliers by shifting.
    # Assume that time is in dimension 1.
    # -----------------------------------------------------------------

    # Need to find the index for each Pg Qg Vm Va using linearindex of the jump model..
    start = linearindex(m_cur[:Pg][1])
    copyvars(inner_cur, start, inner_prev, start+num_gens, num_gens*(T-1))

    start = linearindex(m_cur[:Qg][1])
    copyvars(inner_cur, start, inner_prev, start+num_gens, num_gens*(T-1))

    start = linearindex(m_cur[:Vm][1])
    copyvars(inner_cur, start, inner_prev, start+num_buses, num_buses*(T-1))

    start = linearindex(m_cur[:Va][1])
    copyvars(inner_cur, start, inner_prev, start+num_buses, num_buses*(T-1))

    if piecewise
        start = linearindex(m_cur[:Cg][1])
        copyvars(inner_cur, start, inner_prev, start+num_gens, num_gens*(T-1))
    end

    # -----------------------------------------------------------------
    # Copy constraint values and their multipliers by shifting.
    # Assume that time is in dimension 1.
    # -----------------------------------------------------------------

    if T > 2

        # -------------------------------------------------------------
        # Move the nonzero multipliers of the first ramping constraint
        # to the appropriate bound multipliers of variables.
        #
        # This is because the first ramping constraint was merged into
        # the bound constraints as we shift the time horizon.
        #
        #   e.g., setlower(Pg[2], max(Pg[2].Pmin, Pg[1] - ramping))
        #         setupper(Pg[2], min(Pg[2].Pmax, Pg[1] + ramping))
        #
        # The sign of a multiplier depends on the direction of its
        # constraint:
        #
        #  '>=' has a nonpositive sign.
        #  '<=' has a nonnegative sign.
        #  'a <= g(x) <= b' is equivalent to g(x) >= a and g(x) <= b.
        #   Then the sign follows the previous rule.
        #
        # Also a constraint value depends on the type.
        #
        #   nonlinear constraint: value of g == g(x) - rhs if one-sided,
        #                         value of g == g(x) if double-sided.
        #   linear constraint   : value of g == g(x)
        #
        # For multipliers on variable bounds, they all have nonnegative
        # signs.
        #
        # There is no consistency in JuMP.
        # -------------------------------------------------------------

        Pgstart = linearindex(m_cur[:Pg][1])
        tol = 1e-6

        start = linearindex(m_cur[:ramping][1])
        copyconstrs(inner_cur, start, inner_prev,
                    start+num_gens, num_gens*(T-2))

        # ---------------------------------------------------------
        # Identify active ramping constraints and copy their
        # multipliers to the appropriate bound multipliers.
        # ---------------------------------------------------------

        for g=1:num_gens
            val_g = inner_prev.g[start+g-1]
            mult_g = inner_prev.mult_g[start+g-1]
            ramp_agc = circuit.gen[g].ramp_agc

            if val_g <=  -ramp_agc + tol && mult_g < -tol

                # -------------------------------------------------
                # Ramping down has occured, thus lower bound is
                # active. Copy the multiplier to lower bound mult.
                # -------------------------------------------------

                inner_cur.mult_x_L[Pgstart+g-1] = -mult_g
            elseif val_g >= ramp_agc - tol && mult_g > tol

                # -------------------------------------------------
                # Ramping up has occured, thus upper bound is
                # active. Copy the multiplier to upper bound mult.
                # -------------------------------------------------

                inner_cur.mult_x_U[Pgstart+g-1] = mult_g
            end
        end
    end

    if piecewise
        start = linearindex(m_cur[:plcurve][1,1,1])
        copyconstrs(inner_cur, start, inner_prev, start+size_pl, size_pl*(T-1))
    end

    start = num_linconstr + linearindex(m_cur[:pfreal][1])
    copyconstrs(inner_cur, start, inner_prev, start+num_buses, num_buses*(T-1))

    start = num_linconstr + linearindex(m_cur[:pfimag][1])
    copyconstrs(inner_cur, start, inner_prev, start+num_buses, num_buses*(T-1))

    if num_linelimits > 0
        start = num_linconstr + linearindex(m_cur[:flowmaxfrom][1])
        copyconstrs(inner_cur, start, inner_prev,
                    start+num_linelimits, num_linelimits*(T-1))

        start = num_linconstr + linearindex(m_cur[:flowmaxto][1])
        copyconstrs(inner_cur, start, inner_prev,
                    start+num_linelimits, num_linelimits*(T-1))
    end

    if warm_type == :shift_copy

        # -----------------------------------------------------------------
        # Copy the (T-1)th values for the incoming period.
        # -----------------------------------------------------------------

        start = linearindex(m_cur[:Pg][T,1])
        copyvars(inner_cur, start, inner_prev, start, num_gens)

        start = linearindex(m_cur[:Qg][T,1])
        copyvars(inner_cur, start, inner_prev, start, num_gens)

        start = linearindex(m_cur[:Vm][T,1])
        copyvars(inner_cur, start, inner_prev, start, num_buses)

        start = linearindex(m_cur[:Va][T,1])
        copyvars(inner_cur, start, inner_prev, start, num_buses)

        if piecewise
            start = linearindex(m_cur[:Cg][T,1])
            copyvars(inner_cur, start, inner_prev, start, num_gens)

            start = linearindex(m_cur[:plcurve][T,1,1])
            copyconstrs(inner_cur, start, inner_prev, start, size_pl)
        end

        start = linearindex(m_cur[:ramping][T-1,1])
        copyconstrs(inner_cur, start, inner_prev, start, num_gens)

        # -----------------------------------------------------------------
        # Compute the initial KKT error.
        # -----------------------------------------------------------------

        d_err = max(maximum(abs.(demand.pd[:,T] .- demand.pd[:,T-1])),
                    maximum(abs.(demand.qd[:,T] .- demand.qd[:,T-1])))

        ω_err = 0
        for i=1:num_gens
            if ω_err < abs(inner_prev.mult_g[start + i - 1])
                ω_err = abs(inner_prev.mult_g[start + i - 1])
            end
        end

        println("---------------------------------------------------------")
        println("Initial KKT error: |Δd| = ", d_err, " |ω| = ", ω_err)
        println("---------------------------------------------------------")
        flush(stdout)

        start = num_linconstr + linearindex(m_cur[:pfreal][T,1])
        copyconstrs(inner_cur, start, inner_prev, start, num_buses)

        start = num_linconstr + linearindex(m_cur[:pfimag][T,1])
        copyconstrs(inner_cur, start, inner_prev, start, num_buses)

        if num_linelimits > 0
            start = num_linconstr + linearindex(m_cur[:flowmaxfrom][T,1])
            copyconstrs(inner_cur, start, inner_prev, start, num_linelimits)

            start = num_linconstr + linearindex(m_cur[:flowmaxto][T,1])
            copyconstrs(inner_cur, start, inner_prev, start, num_linelimits)
        end
    elseif warm_type == :shift_phase1

        # -----------------------------------------------------------------
        # Solve a subproblem consisting of only the incoming period.
        # -----------------------------------------------------------------

        p1_demand = Load(demand.pd[:,T], demand.qd[:,T])
        start = linearindex(m_prev[:Pg][T,1])
        prev_val = inner_prev.x[start:start+num_gens-1]
        m_p1 = get_mpcmodel(circuit, p1_demand;
                            has_ramping=true, phase1=true, piecewise=piecewise,
                            prev_val=prev_val)

        # -----------------------------------------------------------------
        # Warm-start the subproblem.
        # -----------------------------------------------------------------

        setsolver(m_p1, IpoptSolver(option_file_name="ipopt.op2"))
        JuMP.build(m_p1)

        in_p1 = internalmodel(m_p1)
        inner_p1 = in_p1.inner
        p1_num_linconstr = MathProgBase.numlinconstr(m_p1)

        start = linearindex(m_cur[:Pg][T,1])
        p1start = linearindex(m_p1[:Pg][1])
        copyvars(inner_p1, p1start, inner_prev, start, num_gens)

        start = linearindex(m_cur[:Qg][T,1])
        p1start = linearindex(m_p1[:Qg][1])
        copyvars(inner_p1, p1start, inner_prev, start, num_gens)

        start = linearindex(m_cur[:Vm][T,1])
        p1start = linearindex(m_p1[:Vm][1])
        copyvars(inner_p1, p1start, inner_prev, start, num_buses)

        start = linearindex(m_cur[:Va][T,1])
        p1start = linearindex(m_p1[:Va][1])
        copyvars(inner_p1, p1start, inner_prev, start, num_buses)

        if piecewise
            start = linearindex(m_cur[:Cg][T,1])
            p1start = linearindex(m_p1[:Cg][1])
            copyvars(inner_p1, p1start, inner_prev, start, num_gens)

            start = linearindex(m_cur[:plcurve][T,1,1])
            p1start = linearindex(m_p1[:plcurve][1,1,1])
            copyconstrs(inner_p1, p1start, inner_prev, start, size_pl)
        end

        start = linearindex(m_cur[:ramping][T-1,1])
        p1start = linearindex(m_p1[:ramping][1])
        copyconstrs(inner_p1, p1start, inner_prev, start, num_gens)

        start = num_linconstr + linearindex(m_cur[:pfreal][T,1])
        p1start = p1_num_linconstr + linearindex(m_p1[:pfreal][1])
        copyconstrs(inner_p1, p1start, inner_prev, start, num_buses)

        start = num_linconstr + linearindex(m_cur[:pfimag][T,1])
        p1start = p1_num_linconstr + linearindex(m_p1[:pfimag][1])
        copyconstrs(inner_p1, p1start, inner_prev, start, num_buses)

        if num_linelimits > 0
            start = num_linconstr + linearindex(m_cur[:flowmaxfrom][T,1])
            p1start = p1_num_linconstr + linearindex(m_p1[:flowmaxfrom][1])
            copyconstrs(inner_p1, p1start, inner_prev, start, num_linelimits)

            start = num_linconstr + linearindex(m_cur[:flowmaxto][T,1])
            p1start = p1_num_linconstr + linearindex(m_p1[:flowmaxto][1])
            copyconstrs(inner_p1, p1start, inner_prev, start, num_linelimits)
        end

        MathProgBase.setwarmstart!(in_p1, inner_p1.x)
        MathProgBase.optimize!(in_p1)

        stat = MathProgBase.status(in_p1)

        if stat != :Infeasible && stat != :Unbounded
            m_p1.colVal = MathProgBase.getsolution(in_p1)
        else
            println("SPOPF has failed with status: ", stat)
            @assert false
        end

        save_rampinfo(m_p1, 2, num_gens, "rampinfo_p1_"*profname; circuit=circuit)

        # -----------------------------------------------------------------
        # Copy the values from the Phase1 solution for the incoming period.
        # -----------------------------------------------------------------

        start = linearindex(m_cur[:Pg][T,1])
        p1start = linearindex(m_p1[:Pg][1])
        copyvars(inner_cur, start, inner_p1, p1start, num_gens)

        start = linearindex(m_cur[:Qg][T,1])
        p1start = linearindex(m_p1[:Qg][1])
        copyvars(inner_cur, start, inner_p1, p1start, num_gens)

        start = linearindex(m_cur[:Vm][T,1])
        p1start = linearindex(m_p1[:Vm][1])
        copyvars(inner_cur, start, inner_p1, p1start, num_buses)

        start = linearindex(m_cur[:Va][T,1])
        p1start = linearindex(m_p1[:Va][1])
        copyvars(inner_cur, start, inner_p1, p1start, num_buses)

        if piecewise
            start = linearindex(m_cur[:Cg][T,1])
            p1start = linearindex(m_p1[:Cg][1])
            copyvars(inner_cur, start, inner_p1, p1start, num_gens)

            start = linearindex(m_cur[:plcurve][T,1,1])
            p1start = linearindex(m_p1[:plcurve][1,1,1])
            copyconstrs(inner_cur, start, inner_p1, p1start, size_pl)
        end

        start = linearindex(m_cur[:ramping][T-1,1])
        p1start = linearindex(m_p1[:ramping][1])
        copyconstrs(inner_cur, start, inner_p1, p1start, num_gens)

        # -----------------------------------------------------------------
        # Compute the initial KKT error.
        # -----------------------------------------------------------------

        ω_err = 0
        for i=1:num_gens
            if ω_err < abs(inner_p1.mult_g[p1start + i - 1])
                ω_err = abs(inner_p1.mult_g[p1start + i - 1])
            end
        end

        println("---------------------------------------------------------")
        println("Initial KKT error: |ω| = ", ω_err)
        println("---------------------------------------------------------")
        flush(stdout)

        start = num_linconstr + linearindex(m_cur[:pfreal][T,1])
        p1start = p1_num_linconstr + linearindex(m_p1[:pfreal][1])
        copyconstrs(inner_cur, start, inner_p1, p1start, num_buses)

        start = num_linconstr + linearindex(m_cur[:pfimag][T,1])
        p1start = p1_num_linconstr + linearindex(m_p1[:pfimag][1])
        copyconstrs(inner_cur, start, inner_p1, p1start, num_buses)

        if num_linelimits > 0
            start = num_linconstr + linearindex(m_cur[:flowmaxfrom][T,1])
            p1start = p1_num_linconstr + linearindex(m_p1[:flowmaxfrom][1])
            copyconstrs(inner_cur, start, inner_p1, p1start, num_linelimits)

            start = num_linconstr + linearindex(m_cur[:flowmaxto][T,1])
            p1start = p1_num_linconstr + linearindex(m_p1[:flowmaxto][1])
            copyconstrs(inner_cur, start, inner_p1, p1start, num_linelimits)
        end
    end
end

function save_header(m, circuit, T, basename="solution_header")
    num_linconstr = MathProgBase.numlinconstr(m)
    num_buses = length(circuit.bus)
    num_gens = length(circuit.gen)
    rateA = getfield.(circuit.line, :rateA)  # getatttr version of struct
    num_linelimits = length(findall((rateA .!= 0) .& (rateA .< 1.0e10)))

    # ----------------------------------------------------------------
    # Save the relative locations.
    # ----------------------------------------------------------------

    f = open(basename*".txt", "w")
    write(f, string(num_linconstr), "\n")
    write(f, string(T), "\n")
    write(f, string(num_gens), "\n")
    write(f, string(num_buses), "\n")
    write(f, string(num_linelimits), "\n")

    start = linearindex(m[:Pg][1])
    write(f, string(start), "\n")

    start = linearindex(m[:Qg][1])
    write(f, string(start), "\n")

    start = linearindex(m[:Vm][1])
    write(f, string(start), "\n")

    start = linearindex(m[:Va][1])
    write(f, string(start), "\n")

    if T > 1
        start = linearindex(m[:ramping][1])
        write(f, string(start), "\n")
    end

    start = num_linconstr + linearindex(m[:pfreal][1])
    write(f, string(start), "\n")

    start = num_linconstr + linearindex(m[:pfimag][1])
    write(f, string(start), "\n")

    if num_linelimits > 0
	start = num_linconstr + linearindex(m[:flowmaxfrom][1])
	write(f, string(start), "\n")

	start = num_linconstr + linearindex(m[:flowmaxto][1])
	write(f, string(start), "\n")
    end

    close(f)
end

function save_rampagc(circuit, basename="solution_rampagc")
    gen = circuit.gen

    f = open(basename*".txt", "w")
    for g = 1:length(gen)
        write(f, string(gen[g].ramp_agc), "\n")
    end
    close(f)
end

function load_solution(m, basename="solution"; noreplace=true)
    if !noreplace
        basename = replace(basename, "cold" => "warm")
        basename = replace(basename, "phase1" => "phase0")
        basename = replace(basename, "use_qp" => "no_qp")
        basename = replace(basename, r"_pf[\d]" => "")
        basename = replace(basename, r"pt0.[\d]+" => "pt0")
        basename = replace(basename, r"_mlim_[\d]+" => "")
    end

    println("Loading a solution from the file: ", basename*".jld2")
    sol = load(basename*".jld2")

    if !m.internalModelLoaded
        JuMP.build(m)
    end

    @assert(sol["n"] == MathProgBase.numvar(m))
    @assert(sol["m"] == MathProgBase.numconstr(m))

    inner = internalmodel(m).inner

    copyto!(inner.x, 1, sol["x"], 1, sol["n"])
    copyto!(inner.mult_x_L, 1, sol["mult_x_L"], 1, sol["n"])
    copyto!(inner.mult_x_U, 1, sol["mult_x_U"], 1, sol["n"])
    copyto!(inner.g, 1, sol["g"], 1, sol["m"])
    copyto!(inner.mult_g, 1, sol["mult_g"], 1, sol["m"])
    m.colVal = MathProgBase.getsolution(internalmodel(m))
    m.objVal = MathProgBase.getobjval(internalmodel(m))
end

function save_solution(m, basename="solution"; dir_name="./results/")
    inner = internalmodel(m).inner
    basename = dir_name * basename

    save(basename*".jld2", "n", length(inner.x), "m", length(inner.g), "x", inner.x,
         "mult_x_L", inner.mult_x_L, "mult_x_U", inner.mult_x_U,
         "g", inner.g, "mult_g", inner.mult_g)
end

function save_rampinfo(m, T, num_gens, basename="rampinfo", mode="a";
		       circuit = nothing, dir_name="./results/")
    if T <= 1
        return
    end

    basename = dir_name * basename
    f = open(basename*".txt", mode)
    inner = internalmodel(m).inner

    # ----------------------------------------------------------------
    # Save the ratio of binded ramping constraints.
    # ----------------------------------------------------------------

    start = linearindex(m[:ramping][1])

    for t=1:T-1
        idx = start + (t-1)*num_gens
        num_binded = 0
        num_constr_binded = 0

        for g=1:num_gens
            if abs(inner.mult_g[idx + g - 1]) >= 1e-3
                num_binded += 1  # where ramping is actually larger than 0
            end

            if circuit != nothing
                if abs(abs(inner.g[idx + g - 1]) - circuit.gen[g].ramp_agc) < 1e-6
                    num_constr_binded += 1  # Binded to the bounds
                end
            end
        end

        s = @sprintf("(%.6f  %.6f)", num_binded/num_gens, num_constr_binded/num_gens)

        if t == 1
            write(f, s)
        else
            write(f, "\t", s)
        end
    end

    write(f, "\n")
    close(f)
end

function save_info(info_dict, profname; dir_name="./results/")
    fname = profname * ".jld2"
    fname = dir_name * fname
    FileIO.save(fname, info_dict)

    println("Arguments are being saved in $(fname)")
    return
end

function solve_cold(m, circuit, demand; powerflow_solve = false)
    max_iter = 1000
    stat = :Infeasible

    if powerflow_solve
        m_pf = get_mpcpfmodel(circuit, demand)
        setsolver(m_pf, IpoptSolver(option_file_name="ipopt.opt"))
        stat = solve(m_pf)
        if stat != :Optimal
            println("Power flow stat is not optimal: ", stat)
            return
        end

        setsolver(m, IpoptSolver(option_file_name="ipopt.pf"))
        JuMP.build(m)
        MathProgBase.setwarmstart!(internalmodel(m), m_pf.colVal)
        MathProgBase.optimize!(internalmodel(m))
        stat = MathProgBase.status(internalmodel(m))
        if stat != :Infeasible && stat != :Unbounded
            T = size(demand.pd, 2)
            m.colVal = MathProgBase.getsolution(internalmodel(m))
            m.objVal = get_mpcobjectivevalue(m, circuit, T)
        end
    else
        init_x(m, circuit, demand)
        setsolver(m, IpoptSolver(option_file_name="ipopt.opt",
                                 max_iter=max_iter,
#                                 linear_solver="ma27",
                                 linear_scaling_on_demand="no"
                                )
        )
#        JuMP.set_optimizer_attribute(m, "linear_solver", "ma27")
        stat = solve(m)
    end

    return stat
end

function print_args(arg_dict)
    println("---------------Parsed Options-----------------")
    for (k, v) in arg_dict
        println("$(@sprintf("%10s", k)) - $(v)")
    end
    println("----------------------------------------------")
    return
end

function usage()
    println("Usage: julia mpc.jl case scen T H LS RS warm opt [profname]")
    println(" where")
    println("          case - the name of the case file")
    println("          scen - the name of the scenario file")
    println("             T - the length of the time horizon")
    println("             H - the length of the time horizon shift")
    println("            LS - load scale")
    println("            RS - ramping scale")
    println("          warm - cold|shift_copy|shift_phase1")
    println("           opt - option file number for warm-start")
    println("      cut_line - 1 to cut a line, 0 otherwise")
    println("       cut_gen - 1 to turn off a generator, 0 otherwise")
    println("       perturb - perturbation scale")
    println("            qp - 1 if qp approx is used 0 otherwise")
    println("      load_sol - 1 if solution for time horizon 1 is available")
    println("      pf_solve - 1 if want to solve power flow for init point")
    println("      profname - profile name")
end

function main(args)
    # ---------------------------------------------------------------------
    # Parse the arguments.
    # ---------------------------------------------------------------------

    if length(args) < 14
        usage()
        return
    end

    if typeof(args) == Array{Any, 1}
        case = args[1]
        scen = args[2]
        T = max(parse(Int,args[3]),1)
        H = max(parse(Int,args[4]),1)
        load_scale = parse(Float64,args[5])
        ramp_scale = parse(Float64,args[6])
        warm = args[7]
        opt = max(parse(Int,args[8]),2)
        cut_line = (parse(Int,args[9]) == 1) ? true : false
        cut_gen = (parse(Int,args[10]) == 1) ? true : false
        perturb = parse(Float64,args[11])
        qp = parse(Int,args[12])
        load_sol = parse(Int,args[13])
        powerflow_solve = (parse(Int,args[14]) == 1) ? true : false
    elseif typeof(args) == OrderedDict{String, Any}
        case = args["case"]
        scen = args["scen"]
        T = args["T"]
        H = args["H"]
        load_scale = args["load_scale"]
        ramp_scale = args["ramp_scale"]
        warm = args["warm"]
        opt = args["opt"]
        cut_line = args["cut_line"]
        cut_gen = args["cut_gen"]
        perturb = args["perturb"]
        qp = args["qp"]
        load_sol = args["load_sol"]
        powerflow_solve = args["powerflow_solve"]
        cut_time = args["cut_time"]
    end

    profname = nothing
    save_dir = "./results_test/$(basename(case))_T-$(T)_cut-$(cut_line)/"  # FIXME
    if ~ispath(save_dir)
        mkpath(save_dir)
    end

    warmstart = true
    if warm == "cold"
        warmstart = false
    elseif warm == "warm"
        warmstart = true
    elseif warm == "shift_copy"
        warm_type = :shift_copy
    elseif warm == "shift_phase1"
        warm_type = :shift_phase1
    else
        usage()
        return
    end

    if cut_line == 1 && cut_gen == 1
        println("Error: a line and a generator cannot be both turned off.")
        usage()
        return
    end

    baseMVA = 100
    piecewise_cost = false

    if profname == nothing
        profname = get_profname(case, T, load_scale, ramp_scale, warm)
    end

    # TODO Redundant code. Need refactoring
    arg_dict = OrderedDict()
    arg_dict["case"] = case
    arg_dict["scen"] = scen
    arg_dict["T"] = T
    arg_dict["H"] = H
    arg_dict["load_scale"] = load_scale  # Can multiply by uniform 0.X ~ 1.X but need range too
    arg_dict["ramp_scale"] = ramp_scale  # ramp_agc = gen.Pmax * ramp_scaling
    arg_dict["warm"] = warm
    arg_dict["opt"] = opt
    arg_dict["cut_line"] = cut_line
    arg_dict["cut_gen"] = cut_gen
    arg_dict["perturb"] = perturb
    arg_dict["qp"] = qp
    arg_dict["load_sol"] = load_sol
    arg_dict["pf_solve"] = powerflow_solve
    arg_dict["profname"] = profname
    arg_dict["cut_time"] = cut_time

    print_args(arg_dict)

    info_dict = OrderedDict()
    info_dict["args"] = arg_dict
    info_dict["cut_time"] = cut_time

    # ---------------------------------------------------------------------
    # Read the circuit and load.
    # circuit, load are instances of Circuit, Load after preprocessing
    # ---------------------------------------------------------------------

    circuit = getcircuit(case, baseMVA, ramp_scale)
    load = getload(scen, circuit, load_scale)

    println("Network statistics:")
    println("   # buses     : ", length(circuit.bus))
    println("   # generators: ", length(circuit.gen))
    println("   # branches  : ", length(circuit.line))

    if length(circuit.gen) > 0 && circuit.gen[1].gentype == 1
        piecewise_cost = true
    end

    println("Number of lines: ", length(circuit.line))
    flush(stdout)

    # ---------------------------------------------------------------------
    # If cut_line is set, we randomly choose from one of the lines
    # ---------------------------------------------------------------------

    cut_circuit = nothing
    if cut_line || cut_gen
        neg_line = neg_gen = -1
#        circuit = getcircuit(case, baseMVA, ramp_scale;
#                             neg_line=neg_line, neg_gen=neg_gen)
        cut_circuit = getcircuit(case, baseMVA, ramp_scale;
                             neg_line="random", neg_gen=neg_gen)
        info_dict["cut_line"] = cut_circuit.cut_line
    end

    println("Number of lines: ", length(circuit.line))
    if cut_line
        println("Number of lines after cut: ", length(cut_circuit.line))
    end
    flush(stdout)

    num_gens = length(circuit.gen)
    gen = circuit.gen

    # ---------------------------------------------------------------------
    # Solve the first time horizon [1:T] using Ipopt.
    # ---------------------------------------------------------------------

    single_objval = 0
    total_objval = 0
    total_single_objval = 0

    demand = Load(load.pd[:,1:T], load.qd[:,1:T])
    m_cur = nothing
    if cut_line  # FIXME version 1, 2
        m_cur = get_cut_mpcmodel2(circuit, cut_circuit, demand, cut_time; piecewise=piecewise_cost)
    else
        m_cur = get_mpcmodel(circuit, demand; piecewise=piecewise_cost)
    end
#    init_x(m_cur, circuit, demand)
#    setsolver(m_cur, IpoptSolver(option_file_name="ipopt.opt"))

    info_dict["Pd"] = demand.pd'  # From 1:T  [1:G, 1:T] -> [1:T, 1:G]
    info_dict["Qd"] = demand.qd'

    if load_sol == 1
        stat = :Optimal
#        load_solution(m_cur, "solution_"*profname*"_h1")
        load_solution(m_cur, "solution_"*profname; noreplace=true)
        save_rampinfo(m_cur, T, num_gens, "rampinfo_"*profname, "w"; circuit=circuit)
    else
        if !warmstart
            stat = solve_cold(m_cur, circuit, demand; powerflow_solve = powerflow_solve)
        else
            if cut_line
                println("Generating warm start solution...........")
                demand0 = Load(load.pd[:,1:cut_time], load.qd[:,1:cut_time])
                m_warm = get_mpcmodel(circuit, demand0; piecewise=piecewise_cost)
                stat = solve_cold(m_warm, circuit, demand0; powerflow_solve = powerflow_solve)
#                m_warm = get_mpcmodel(circuit, demand; piecewise=piecewise_cost)
#                stat = solve_cold(m_warm, circuit, demand; powerflow_solve = powerflow_solve)

#                init_x(m_warm, circuit, demand)
#                setsolver(m_warm, IpoptSolver(option_file_name="ipopt.opt"))
#                stat = solve(m_warm)
                if stat != :Optimal
                    println("Warm start stat is not optimal: ", stat)
                    return
                end

#                JuMP.build(m_cur)
                init_x(m_cur, m_warm, circuit, cut_circuit, demand, cut_time, T)
                fix_vars_until_cut(m_cur, circuit, cut_time)  # FIXME
#                real_viol, imag_viol = get_powerflow_violation(m_cur, circuit, demand; periods=collect(1:cut_time))
#                println(max(vec(real_viol)...))
#                println(max(vec(imag_viol)...))


    #            println(getvalue(m_cur[:Pg])[1:10])
    #            println(getvalue(m_cur[:Pg])[1:10])
            else
                init_x(m_cur, circuit, demand)
            end

            println("Now solving with warm start...")
            setsolver(m_cur, IpoptSolver(option_file_name="ipopt.opt",
                                     max_iter=1000,
#                                     linear_solver="ma27",
                                     linear_scaling_on_demand="no"
                                    )
            )
            stat = solve(m_cur)
        end
        if stat != :Optimal
            println("Stat is not optimal: ", stat)
            return
        end
#        if stat == :Optimal
#            println("Stat is optimal: ", stat, " Passing this example")
#            return
#        else
#            println("Stat is not optimal: ", stat)
#        end

        #### save info......................!
        Pg = getvalue(m_cur[:Pg])
        info_dict["Pg"] = Pg
        Sg = zeros(size(Pg, 1)-1, size(Pg, 2))
        for i in 1:(size(Pg, 1)-1)
            Sg[i, :] = Pg[i+1, :] - Pg[i, :]
        end

        info_dict["Sg"] = Sg
        info_dict["Qg"] = getvalue(m_cur[:Qg])
        info_dict["Vm"] = getvalue(m_cur[:Vm])
        info_dict["Va"] = getvalue(m_cur[:Va])
        save_solution(m_cur, "solution_"*profname; dir_name=save_dir)
        save_rampinfo(m_cur, T, num_gens, "rampinfo_"*profname, "w"; circuit=circuit, dir_name=save_dir)

    end

    single_objval = get_mpcobjectivevalue(m_cur, circuit, 1)
    total_single_objval += single_objval
    total_objval += single_objval

    # ---------------------------------------------------------------------
    # Warm-start option for Ipopt.
    # ---------------------------------------------------------------------
    m_prev = nothing
#    Random.seed!(0)

    # Add the remaining cost.
    total_objval += getobjectivevalue(m_cur) - single_objval

    @printf("Total objective value..............: %18.6e\n", total_objval)
    @printf("Total single objective value.......: %18.6e\n", total_single_objval)

    # Saving all information

    in_cur = internalmodel(m_cur)
    info_dict["total_objval"] = total_objval
    info_dict["total_single_objval"] = total_single_objval
    info_dict["p_status"] = MathProgBase.status(in_cur)
    save_info(info_dict, "info_"*profname; dir_name=save_dir)
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    println(ARGS)
    println("INFO: SET h=1 for single horizon MPOPF")

    if length(ARGS) == 15  # Not a good way to make arguments. Try using ArgParse or use a Dict
        niter = max(parse(Int,ARGS[15]),1)
    end

    for i in 1:niter
        main(ARGS)
    end
end


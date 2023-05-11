using Printf
using Distributions

mutable struct Bus
    # each column in the case**.bus matches these values in sequential order.
    # each row represents a different bus
    bus_i::Int
    bustype::Int
    Pd::Float32
    Qd::Float32
    Gs::Float32
    Bs::Float32
    area::Int
    Vm::Float32
    Va::Float32
    baseKV::Float32
    zone::Int
    Vmax::Float32
    Vmin::Float32
end

mutable struct Line
    # impedance z = r + jx
    # admittance y = g + jb
    from::Int
    to::Int
    r::Float32     # r: resistance
    x::Float32     # x: reactance
    b::Float32     # b: susceptance
    rateA::Float32
    rateB::Float32
    rateC::Float32
    ratio::Float32 #TAP
    angle::Float32 #SHIFT
    status::Int
    angmin::Float32
    angmax::Float32
end

mutable struct Gen
    # .gen fields
    bus::Int
    Pg::Float32
    Qg::Float32
    Qmax::Float32
    Qmin::Float32
    Vg::Float32
    mBase::Float32
    status::Int
    Pmax::Float32
    Pmin::Float32
    Pc1::Float32
    Pc2::Float32
    Qc1min::Float32
    Qc1max::Float32
    Qc2min::Float32
    Qc2max::Float32
    ramp_agc::Float32      # Equivalent to r_g
    # .gencost fields
    gentype::Int
    startup::Float32
    shutdown::Float32
    n::Int
    coeff::Array{Float32}  # Cost matrix for a generator that is turned "on"
end

struct Yline
    # Impedance matrix
    # Y = / YffR + j YffI, YftR + j YftI \
    #     \ YtfR + j YtfI, YttR + j YttI /
    from::Int
    to::Int
    YffR::Float32  # g^c_{ik}
    YffI::Float32  # b^c_{ik}
    YttR::Float32  # g^c_{ki}
    YttI::Float32  # b^c_{ki}
    YtfR::Float32  # g_{ki}
    YtfI::Float32  # b_{ki}
    YftR::Float32  # g_{ik}
    YftI::Float32  # b_{ik}
end

struct Ybus
    bus::Int
    YshR::Float32  # g^s_i
    YshI::Float32  # b^s_i
end

mutable struct Circuit
    baseMVA::Float32
    busref::Int
    bus::Array{Bus}
    line::Array{Line}
    gen::Array{Gen}
    yline::Array{Yline}
    ybus::Array{Ybus}
    busdict::Dict{Int,Int}
    frombus::Array
    tobus::Array
    bus2gen::Array
    cut_line::Int
    is_bus2gen::Array{Float32, 2}  # is_bus2gen[b, i] = 1 if bus i \in bus2gen[b] else 0
    bar_YftR::Array{Float32, 2}  # (b, i) for i in 1:n_bus
    bar_YftI::Array{Float32, 2}
    bar_YtfR::Array{Float32, 2}  # (i, b)
    bar_YtfI::Array{Float32, 2}
    bar_YshR::Array{Float32, 1}  # (b)
    bar_YshI::Array{Float32, 1}  # (b)

    Circuit(baseMVA, busref, bus, line, gen, yline, ybus, busdict,
            frombus, tobus, bus2gen, cut_line) = new(
        baseMVA, busref, bus, line, gen, yline, ybus, busdict,
        frombus, tobus, bus2gen, cut_line,
        get_is_bus2gen(bus2gen, bus, gen),
        get_bar_Y(bus, ybus, frombus, tobus, yline, line)...
    )
end

struct Load
    # pd - demand for p (real)
    # qd - demand for q (imaginary)
    pd
    qd
end

function get_is_bus2gen(bus2gen, bus, gen)
    n_bus = length(bus)
    n_gen = length(gen)
    return [Float32(g in bus2gen[b]) for b=1:n_bus, g=1:n_gen]
#    is_bus2gen = zeros(Float32, (n_bus, n_bus))
#    for b in 1:n_bus
#        for i in 1:n_bus
#            is_bus2gen[b, i] = i in bus2gen[b] ? 1 : 0
#        end
#    end
#
end

function get_bar_Y(bus, ybus, frombus, tobus, yline, line)
    n_bus = length(bus)
    bar_YftR = zeros(Float32, (n_bus, n_bus))  # bar_YftR[b, i]
    bar_YftI = zeros(Float32, (n_bus, n_bus))
    bar_YtfR = zeros(Float32, (n_bus, n_bus))  # bar_YtfR[i, b]
    bar_YtfI = zeros(Float32, (n_bus, n_bus))

    for b in 1:n_bus
        line_to = [line[l].to for l in frombus[b]]
        for i in 1:n_bus
            if i in line_to
                l_idx = findall(x -> x == i, line_to)[1]
                l = frombus[b][l_idx]
                bar_YftR[b, i] = yline[l].YftR
                bar_YftI[b, i] = yline[l].YftI
            else
                bar_YftR[b, i] = 0
                bar_YftI[b, i] = 0
            end
        end
    end

    for b in 1:n_bus
        line_from = [line[l].from for l in tobus[b]]
        for i in 1:n_bus
            if i in line_from
                l_idx = findall(x -> x == i, line_from)[1]
                l = tobus[b][l_idx]
                bar_YtfR[i, b] = yline[l].YtfR
                bar_YtfI[i, b] = yline[l].YtfI
            else
                bar_YtfR[i, b] = 0
                bar_YtfI[i, b] = 0
            end
        end
    end

    bar_YshR = Float32.([
        (sum([yline[l].YffR for l in frombus[b]])
         + sum([yline[l].YttR for l in tobus[b]])
         + ybus[b].YshR)
         for b in 1:n_bus
    ])

    bar_YshI = Float32.([
         (sum([-yline[l].YffI for l in frombus[b]])
         + sum([-yline[l].YttI for l in tobus[b]])
         - ybus[b].YshI)
         for b in 1:n_bus
    ])

    return bar_YftR, bar_YftI, bar_YtfR, bar_YtfI, bar_YshR, bar_YshI
end

function get_busmap(bus)
    # Dict{K, V}
    # Hash table with key type K and value of type V
    busdict = Dict{Int,Int}()

    for i in 1:length(bus)
        # check if busdict has no key with value of bus_i
        @assert !haskey(busdict,bus[i].bus_i)
        busdict[bus[i].bus_i] = i
    end

    return busdict
end

function get_linetobusmap(bus, line, busdict)
    num_buses = length(bus)
    # Array of Arrays Int[]
    from = [Int[] for i in 1:num_buses]
    to = [Int[] for i in 1:num_buses]

    for i in 1:length(line)
        idx = busdict[line[i].from]  # idx of bus that line i comes from
        @assert 1 <= idx <= num_buses
        push!(from[idx], i)

        idx = busdict[line[i].to]
        @assert 1 <= idx <= num_buses
        push!(to[idx], i)
    end

    # from[1] = 2 means that line[2] begins FROM bus 1 where bus 1 is just a number assignment of buses
    # to[3] = 10 means that line[10] goes TO bus 3
    # So it is actually bus to line map..
    # Given a bus B, from[B] is the line number (idx) that starts at B
    #                to[B]   is the line number (idx) that ends a
    #
    #  --------> [bus B] --------->
    #    to[B]            from[B]
    return from, to
end

function get_bustogenmap(bus, gen, busdict)
    bus2gen = [Int[] for i in 1:length(bus)]

    for i in 1:length(gen)
        idx = busdict[gen[i].bus]
        push!(bus2gen[idx], i)
    end

    # Given a bus idx B, bus2gen[B] is the gen number (idx) that contains bus B
    return bus2gen
end

# -------------------------------------------------------------------------
# Compute admittances.
# -------------------------------------------------------------------------
function getY(case, line, bus, baseMVA)
    dim1 = size(line,1)
    Ys = complex(zeros(dim1, 1))
    tap = complex(ones(dim1, 1))
    Ytt = complex(zeros(dim1, 1))
    Yff = complex(zeros(dim1, 1))
    Yft = complex(zeros(dim1, 1))
    Ytf = complex(zeros(dim1, 1))

    # ---------------------------------------------------------------------
    # bus f: tap bus, bus t: impedance bus or Z bus
    #
    # Ys: the admittance between bus f and bus t.
    #     It is the reciprocal of a series impedance between bus f and t.
    #
    # When there is a off-nominal transformer, the admittance matrix is
    # defined as follows:
    #
    #   / Ift \ = / Yff  Yft \ / Vf \
    #   \ Itf / = \ Ytf  Ytt / \ Vt /
    #
    # where
    # (
    #  see skm.com/faq_ptw16.html,
    #      https://www.gridpack.org/wiki/images/7/7e/Ybus.pdf
    #  )
    #
    #    Yff = ((Ys + j*bft/2) / |a|^2)
    #    Yft = (-Ys / conj(a))
    #    Ytf = (-Ys / a)
    #    Ytt = (Ys + j*bft/2)
    #
    # When we compute If or It (total current injection at bus f or t),
    # we need to add the bus shunt, YshR and YshI.
    # ---------------------------------------------------------------------

    Ys = 1 ./ (line[:,3] .+ line[:,4] .* im)  # r  line[:3] = r, line[:4] = x
    ind = findall(line[:,9] .!= 0)
    tap[ind] = line[ind,9]                    # ratio
    tap .*= exp.(line[:,10] .* pi/180 .* im)  # angle
    Ytt = Ys .+ line[:,5] ./ 2 .* im          # b
    Yff = Ytt ./ (tap.*conj.(tap))
    Yft = -Ys ./ conj.(tap)
    Ytf = -Ys ./ tap

    yline = zeros(dim1, 10)
    yline[:,1:2] = line[:,1:2]
    yline[:,3] = real.(Yff); yline[:,4] = imag.(Yff)
    yline[:,5] = real.(Ytt); yline[:,6] = imag.(Ytt)
    yline[:,7] = real.(Ytf); yline[:,8] = imag.(Ytf)
    yline[:,9] = real.(Yft); yline[:,10] = imag.(Yft)

    ybus = zeros(size(bus,1), 3)
    YshR = bus[:,5] ./ baseMVA      # Gs
    YshI = bus[:,6] ./ baseMVA      # Bs
    ybus = [ bus[:,1] YshR YshI ]

    @assert 0==length(findall(isnan.(yline[:,3:10])))
    @assert 0==length(findall(isinf.(yline[:,3:10])))
    @assert 0==length(findall(isnan.(ybus[:,2:3])))
    @assert 0==length(findall(isinf.(ybus[:,2:3])))

    nylines = size(yline,1)
    nybuses = size(ybus,1)
    Ylines = Array{Yline}(undef, nylines)
    Ybuses = Array{Ybus}(undef, nybuses)

    for i in 1:nylines
        Ylines[i] = Yline(yline[i,1:end]...)
    end

    for i in 1:nybuses
        Ybuses[i] = Ybus(ybus[i,1:end]...)
    end

    #TODO: What are Ylines and Y buses?
    return Ylines, Ybuses
end

# -------------------------------------------------------------------------
# Get circuit and power demand data for OPF computation.
# -------------------------------------------------------------------------
function getcircuit(case, baseMVA, ramp_scaling; neg_line=-1, neg_gen=-1)
    #=
    Args:
        case -  name of case in data folder
        baseMVA - TODO
    =#

    # DelimitedFiles.readdlm (read dlm)
    # Reads a matrix from the source where each line (separated by eol) gives one row.
    bus_mat = readdlm(case*".bus")
    bus_mat[:,9] *= pi/180  # multiply Va by pi/180.

    branch_mat = readdlm(case*".branch")
    active_line_ind = findall(branch_mat[:,11] .> 0)

    if neg_line == "random"
        # FIXME need to randomly cut line with probability p and set neg_line = -1
        neg_line = rand(1:length(active_line_ind))
    end

    if neg_line != -1
        l = active_line_ind[neg_line]
        @printf("The line (%5d,%5d) has been cut.\n",
                 branch_mat[l,1], branch_mat[l,2])
        println("Line at index $(l) has been cut")
        # remove item at given index [neg_line] (cut line)
        deleteat!(active_line_ind, neg_line)
    end

    line_mat = branch_mat[active_line_ind,:]

    gen_mat = readdlm(case*".gen")
    gencost_mat = readdlm(case*".gencost")

    # Allocate space for coefficients of a quadratic objective function.
    gens_on = findall(x-> x != 0, gen_mat[:,8])
    if neg_gen != -1
        @printf("Generator at %d connected to bus %d has been turned off.\n",
                 gens_on[neg_gen], gen_mat[gens_on[neg_gen],1])
        deleteat!(gens_on, neg_gen)
    end
    num_on = length(gens_on)

    # Compute admittances.
    yline, ybus = getY(case, line_mat, bus_mat, baseMVA)

    num_buses = size(bus_mat,1)
    num_lines = size(line_mat,1)
    num_gens = size(gen_mat,1)

    # Array of mutable struct
    bus = Array{Bus}(undef, num_buses)  # dim - num_uses x 1
    line = Array{Line}(undef, num_lines)
    gen = Array{Gen}(undef, num_on)

    # assigning elements and check # of reference bus
    busref = -1
    for i in 1:num_buses
        # splat operator equivalent to *bus_mat[i, 1:end] returning a tuple ()
        # or **kwargs for Dict object
        # https://stackoverflow.com/questions/33836360/what-do-the-triple-dots-in-a-julia-array-do-why-do-they-change-the-type-s
        bus[i] = Bus(bus_mat[i,1:end]...)
        if bus[i].bustype == 3
            if busref > 0
                error("More than one reference bus present")
            else
                busref = i
            end
        end
    end

    for i in 1:num_lines
        line[i] = Line(line_mat[i,1:end]...)
    end

    # can use enumerate(gens_on) instead
    # Assigning values manually; division of baseMVA in denominator.
    j = 0
    for i in gens_on
        j += 1

        gen[j] = Gen(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Array{Int}(undef,0))
        gen[j].bus = gen_mat[i,1]
        gen[j].Pg = gen_mat[i,2] / baseMVA
        gen[j].Qg = gen_mat[i,3] / baseMVA
        gen[j].Qmax = gen_mat[i,4] / baseMVA
        gen[j].Qmin = gen_mat[i,5] / baseMVA
        gen[j].Vg = gen_mat[i,6]
        gen[j].mBase = gen_mat[i,7]
        gen[j].status = gen_mat[i,8]
        @assert gen[j].status == 1
        gen[j].Pmax = gen_mat[i,9] / baseMVA
        gen[j].Pmin = gen_mat[i,10] / baseMVA
        gen[j].Pc1 = gen_mat[i,11]
        gen[j].Pc2 = gen_mat[i,12]
        gen[j].Qc1min = gen_mat[i,13]
        gen[j].Qc1max = gen_mat[i,14]
        gen[j].Qc2min = gen_mat[i,15]
        gen[j].Qc2max = gen_mat[i,16]
        if gencost_mat[i,1] == 1
            gen[j].ramp_agc = gen_mat[i,17] / baseMVA
        else
            gen[j].ramp_agc = gen[j].Pmax * ramp_scaling
        end
        gen[j].gentype = gencost_mat[i,1]
        gen[j].startup = gencost_mat[i,2]
        gen[j].shutdown = gencost_mat[i,3]
        gen[j].n = gencost_mat[i,4]
        gen[j].coeff = gencost_mat[i,5:end]
    end

    # Create dictionaries due to the lack of set data structure in JuMP.
    busdict = get_busmap(bus) # busdict[bus_i] = some index between 1:nbus
    frombus, tobus = get_linetobusmap(bus, line, busdict)
    bus2gen = get_bustogenmap(bus, gen, busdict)

    # Circuit instance contrains all the information
    circuit = Circuit(baseMVA, busref, bus, line, gen, yline, ybus,
                      busdict, frombus, tobus, bus2gen, neg_line)

    return circuit
end

function getload(scen, circuit, load_scale=0)  # sample from Uniform distribution with size(pd_mat)...
    pd_mat = readdlm(scen*".Pd")
    qd_mat = readdlm(scen*".Qd")

    if load_scale >= 1 || load_scale < 0
        error("ERROR: invalid load_scale: 0 <= $(load_scale) < 1")
    end

    if load_scale > 0
        pd_mat .*= rand(Uniform(1 - load_scale, 1 + load_scale), size(pd_mat))
        qd_mat .*= rand(Uniform(1 - load_scale, 1 + load_scale), size(qd_mat))
    end

    return Load(pd_mat, qd_mat)
end

using PowerModels
using Random

function write_file(arr, dir)
    io = open(dir, "w") do io
        for x in arr
            println(io, x)
        end
    end
end

# This generates random minimum uptime and downtime and startup/shutdown costs, assuming all units are on at the beginning
# Example:
# mfile_dir = "data/case118.m"
# output_dir = "data/multiperiod_data/case118_gen"
function uc_data_generator_initial_all_on(mfile_dir, output_dir; seed = 0)
    Random.seed!(seed)
    data = PowerModels.parse_file(mfile_dir)
    ngen = sum([v["gen_status"] for (k,v) in data["gen"]])
    
    v0 = ones(Int, ngen)
    tu = rand(0:6, ngen)
    td = rand(0:6, ngen)
    hu = Int[rand(0:tu[i]÷2) for i in 1:ngen]
    hd = zeros(Int, ngen)
    con = Float64[ rand() * 3000 + 2000 for i in 1:ngen]
    coff = Float64[ rand() * 2000 + 1000 for i in 1:ngen]
    
    write_file(v0, "$(output_dir).v0")
    write_file(tu, "$(output_dir).tu")
    write_file(td, "$(output_dir).td")
    write_file(hu, "$(output_dir).hu")
    write_file(hd, "$(output_dir).hd")
    write_file(con, "$(output_dir).con")
    write_file(coff, "$(output_dir).coff")    
end

# This generates random minimum uptime and downtime and startup/shutdown costs, assuming all units are on at the beginning
# Example:
# mfile_dir = "data/case118.m"
# output_dir = "data/multiperiod_data/case118_gen"
function uc_data_generator(mfile_dir, output_dir; seed = 0)
    Random.seed!(seed)
    data = PowerModels.parse_file(mfile_dir)
    ngen = sum([v["gen_status"] for (k,v) in data["gen"]])
    
    v0 = rand(0:1, ngen)
    tu = rand(0:6, ngen)
    td = rand(0:6, ngen)
    hu = Int[rand(0:tu[i]÷2) for i in 1:ngen]
    hu = zeros(Int, ngen)
    hd = zeros(Int, ngen)
    for i in 1:ngen
        if v0[i] == 1
            hu[i] = rand(0:tu[i]÷2)
        # else
        #     hd[i] = rand(0:td[i]÷2) # this may lead to infeasible cases
        end
    end
    con = Float64[ rand() * 3000 + 2000 for i in 1:ngen]
    coff = Float64[ rand() * 2000 + 1000 for i in 1:ngen]
    
    write_file(v0, "$(output_dir).v0")
    write_file(tu, "$(output_dir).tu")
    write_file(td, "$(output_dir).td")
    write_file(hu, "$(output_dir).hu")
    write_file(hd, "$(output_dir).hd")
    write_file(con, "$(output_dir).con")
    write_file(coff, "$(output_dir).coff")    
end

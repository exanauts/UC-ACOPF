include("uc_data_generator_utils.jl")

# mfile_dir = "data/case9.m"
# output_dir = "data/multiperiod_data/case9_gen_initial_all_on"
# output_dir = "data/multiperiod_data/case9_gen"

# mfile_dir = "data/case118.m"
# output_dir = "data/multiperiod_data/case118_gen_initial_all_on"
# output_dir = "data/multiperiod_data/case118_gen"

cases = ["case9", "case118", "case2868rte"]

for case in cases
    mfile_dir = "data_v2/$(case).m"
    output_dir = "data_v2/multiperiod_data/$(case)_gen"
    uc_data_generator_v2(mfile_dir, output_dir)
end
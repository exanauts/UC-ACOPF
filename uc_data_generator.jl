include("uc_data_generator_utils.jl")

mfile_dir = "data/case2868rte.m"
output_dir = "data/multiperiod_data/case2868rte_gen"

uc_data_generator_initial_all_on(mfile_dir, output_dir)
# uc_data_generator(mfile_dir, output_dir)
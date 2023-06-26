using DelimitedFiles
filedir = ARGS[1]
ratio = parse(Float64, ARGS[2])
mat = readdlm(filedir)
mat .*= ratio
writedlm(filedir, mat)
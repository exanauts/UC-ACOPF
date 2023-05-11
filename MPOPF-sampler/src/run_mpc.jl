include("mpc.jl")
include("args.jl")


if abspath(PROGRAM_FILE) == @__FILE__
    println(ARGS)  # (case#, numiter)

    case = ARGS[1]
    niter = max(parse(Int, ARGS[2]), 1)  # number of samples to generate

    if case == "case9"
        args = case9_args
    elseif case == "case30"
        args = case30_args
    elseif case == "case118"
        args = case118_args
    elseif case == "case300"
        args = case300_args
    elseif case == "case1354pegase"
        args = case1354pegase_args
    else
        error("Argument Missing: $(case) unknown. Check args.jl")
    end

    for i in 1:niter
        main(args)
    end
end

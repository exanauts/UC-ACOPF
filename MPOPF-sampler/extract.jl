import FileIO

BASE_PATH = abspath(@__DIR__)

function get_case_sample_data(data_dir="sample_results/case9", data_fname="info_ipopt_case9.jld2")
    pth2file = joinpath(BASE_PATH, data_dir, data_fname)
    info_dict = FileIO.load(pth2file)

    return info_dict
end

function get_case_data(raw_data_dir, data_fname)
    # Reads and loads train, test data for supervised learning
    pth2files = joinpath(BASE_PATH, raw_data_dir)
    save_path = joinpath(BASE_PATH, "results/$(data_fname).jld2")

    cut_lines = nothing
    if ispath(save_path)
        println("Loading preprocessed data from $(save_path)........")
        data_dict = FileIO.load(save_path)
        x_data, y_data = data_dict["x_data"], data_dict["y_data"]
        cut_lines = haskey(data_dict, "cut_lines") ? data_dict["cut_lines"] : repeat([-1], size(x_data, 2))
    else
        println("No data in path. Creating data from $(pth2files)........")
        x_data, y_data, cut_lines = make_save_data(pth2files, save_path)
    end

    return x_data, y_data, cut_lines
end

function make_save_data(pth2files, save_path)
    pth2files = joinpath(BASE_PATH, pth2files)
    fnames = readdir(pth2files)
    fnames = filter(x -> startswith(x, "info"), fnames)
    x_data = nothing
    y_data = nothing
    cut_lines = nothing

    println("Loading individual trajectory info from: $(pth2files)...")
    println("Number of expert trajectories saved as data: $(length(fnames))")

    for (i, fname) in enumerate(fnames)
        println("Processing $(i)th trajectory")
        case_dict = FileIO.load("$(pth2files)/$(fname)")
        if String(case_dict["p_status"]) != "Optimal"
            println("Passing $(case_dict["p_status"])")
            continue
        end
        x_, y_, cut_line = extract_data(case_dict)
        x_data = (i == 1) ? x_ : hcat(x_data, x_)
        y_data = (i == 1) ? y_ : hcat(y_data, y_)
        cut_lines = (i == 1) ? cut_line : vcat(cut_lines, cut_line)

        if i > 5000
            break
        end
    end

    println("Saving preprocessed data to $(save_path)........")
    data_dict = Dict("x_data" => x_data, "y_data" => y_data, "cut_lines" => cut_lines)
    FileIO.save(save_path, data_dict)

    return x_data, y_data, cut_lines
end

function extract_data(case_dict)
    # Extracts data from t=1,...,T in each trajectory
    x_data = nothing
    y_data = nothing
    cut_lines = nothing

    for t in 1:(size(case_dict["Pg"], 1)-1)
        x_datum = [                    # e.g. for case 9
            case_dict["Pg"][t, :]...   # (4,)  1
            case_dict["Qg"][t, :]...   # (4,)  5
            case_dict["Vm"][t, :]...   # (9,)  9
            case_dict["Va"][t, :]...   # (9,)  18
            case_dict["Pd"][t+1, :]... # (9,)  27
            case_dict["Qd"][t+1, :]    # (9,)  36 ~ 44
        ]
        y_datum = [
           (case_dict["Pg"][t+1, :] .- case_dict["Pg"][t, :])... # (4,)  1
            case_dict["Qg"][t+1, :]... # (4,)  5
            case_dict["Vm"][t+1, :]... # (9,)  9
            case_dict["Va"][t+1, :]    # (9,)  18 ~ 26
        ]

        x_data = t == 1 ? x_datum : hcat(x_data, x_datum)
        y_data = t == 1 ? y_datum : hcat(y_data, y_datum)
    end

    cut_line = haskey(case_dict, "cut_line") ? case_dict["cut_line"] : -1
    println("Cut line: $(cut_line)")
#    println("Cut time: $(case_dict["cut_time"])")
    println("Problem solve status: $(case_dict["p_status"])")
    cut_lines = repeat([cut_line], size(case_dict["Pg"], 1)-1)

    if haskey(case_dict, "cut_time") && case_dict["cut_time"] > 1  # Full graph until cut_time
        cut_lines[1:max(case_dict["cut_time"]-1, 1)] .= -1
    end

    return x_data, y_data, cut_lines
end


if abspath(PROGRAM_FILE) == @__FILE__
#    info_dict = get_case_sample_data()
#    println(keys(info_dict))
#    get_case_data("results/case9", "case9_aggregated")

#    cut = true
#    for case_num in [9, 30, 118, 300]
#        case_name = "case$(case_num)_T-30_cut-$(cut)"
#        pth2files = "results/$(case_name)"
#        save_path = "results/$(case_name).jld2"
#        make_save_data(pth2files, save_path)
#    end

    case_num = parse(Int, ARGS[1])
    cut = ARGS[2] == "true" ? true : false

    dir_name = "test"
    case_name = "case$(case_num)_T-30_cut-$(cut)"
    pth2files = "results_$(dir_name)/$(case_name)"
    save_path = "results_$(dir_name)/$(case_name).jld2"
    make_save_data(pth2files, save_path)

end

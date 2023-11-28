using HDF5: h5open
using DelimitedFiles: writedlm, readdlm
# using Interpolations, QuadGK
using Romberg 


export get_ghz_value,
       get_T_value,
       get_lmd_value,
       get_kernel_data, 
       calculate_T1,
       process_file


get_period(ghz::AbstractString) = π * 2.0 * 100. / parse(Float64, ghz)
get_period(ghz::Number) = π * 2.0 * 100. / ghz
rel_devi(a, b) = (a - b) / (a + b) * 2.0
get_lmd(lmd_dir) = parse(Float64, split(lmd_dir, "_")[end])
get_extension(filedir) = split(basename(filedir), ".")[end]

function vec2range(x::Vector)
    xrange = range(minimum(x), maximum(x), length=length(x))
    if collect(xrange) ≈ x
        return xrange
    end
    error("When you want to convert a vector to range, you have input a vector is not evenly spaced!")
end

function get_ghz_value(dir)
    # Get the parent directory
    abs_dir = abspath(dir)

    # Check if the directory matches the pattern "*-GHz-equiv"
    mm1 = match(r"\d+-GHz-equiv", abs_dir)
    if mm1 !== nothing
        # Extract the GHz value
        ghz_value = parse(Int, match(r"\d+", mm1.match).match)
        return ghz_value
    else
        error("Directory does not match the pattern '*-GHz-equiv'.")
    end
end

function get_T_value(dir)
    # Get the parent directory
    abs_dir = abspath(dir)
    
    # Check if the directory matches the pattern "T_*"
    mm1 = match(r"T_\d+.\d+", abs_dir)
    if mm1 !== nothing
        # Extract the GHz value
        temp_value = parse(Float64, match(r"\d+.\d+", mm1.match).match)
        return temp_value
    else
        error("Directory does not match the pattern 'T_\\d+.\\d+'.")
    end
end

function get_lmd_value(dir)
    # Get the parent directory
    abs_dir = abspath(dir)
    
    # Check if the directory matches the pattern "*-GHz-equiv"
    mm1 = match(r"lmd_\d+.\d+", abs_dir)
    if mm1 !== nothing
        # Extract the λ value
        λ_value = parse(Float64, match(r"\d+.\d+", mm1.match).match)
        return λ_value
    else
        error("Directory does not match the pattern 'lmd_\\d+.\\d+'.")
    end
end


function get_mask(x_data, y_data, T; early_truncate=false)
    if (T > x_data[end])
        error("Quadratic kernels are oscillatory. Thus you have to provide an origianl data longer than 1 period.")
    else
        dvdval = x_data[end]/T
        fval = floor(Integer, dvdval)
        cval = ceil(Integer, dvdval)
        maskval = isapprox(dvdval, cval; atol=1e-3) ? cval : fval

    	return early_truncate ? (x_data .< 1.0 * T) : (x_data .< maskval*T)
    end
end


function read_h5_kernel(filename)
    x_data = nothing  # Initialize x_data
    y_data = nothing  # Initialize y_data
    h5open(filename) do f
        x_data = read(f["kernel"]["time"]) 
        y_data = read(f["kernel"]["k1_real"])
    end
    return (x_data, y_data)
end

function read_dat_kernel(filename)
    res = readdlm(filename)
    x_data = res[:, 1]
    y_data = res[:, 2]
    return (x_data, y_data)
end

function get_kernel_data(filename, ghz; early_truncate=false, mask=true)
    
    if get_extension(filename) == "h5"
        (x_data, y_data) = read_h5_kernel(filename)
    elseif get_extension(filename) == "dat"
        (x_data, y_data) = read_dat_kernel(filename)
        tol = 1e-6
	    devi = abs(y_data[end] / y_data[1])
	    if (devi < tol)
	        return x_data, y_data
	    else
	        error("'prop-kernel-eq.dat' not converged. Please use dmd to extend; or you can bump up t_final")
	    end
    else
        # println(get_extension(filename))
        error("You data file isn't 'prop-kernel-eq.dat' or 'dmd.h5'. Check your script! Abort!")
    end    

    T = get_period(ghz)
    if mask
    	mask = get_mask(x_data, y_data, T; early_truncate)
    	x_data = x_data[mask]
    	y_data = y_data[mask]
    end

    if x_data === nothing || y_data === nothing
        error("Failed to read data from file.")
    end

    return (x_data, y_data)
end

function calculate_T1(x_data, y_data)
    xrange = vec2range(x_data)
    result, error = romberg(xrange, y_data) 
    # interp = LinearInterpolation(x_data, y_data)
    # result, error = quadgk(interp, minimum(x_data), maximum(x_data))
    println("Numerical integration has error: ", error)
    T1 = 1.0 / abs(result) * 0.5
    println("T1 time is: ", T1, ".")
    return T1
end

function calculate_T1(x_data_dmd, y_data_dmd, x_data_org, y_data_org)
    xrange_org = vec2range(x_data_org)
    result_org, error_org = romberg(xrange_org, y_data_org) 
    # interp_org = LinearInterpolation(x_data_org, y_data_org)
    # result_org, error_org = quadgk(interp_org, minimum(x_data_org), maximum(x_data_org))

    mask_dmd = (x_data_dmd .>= x_data_org[end])
    x_data_dmd = x_data_dmd[mask_dmd]
    y_data_dmd = y_data_dmd[mask_dmd]
    xrange_dmd = vec2range(x_data_dmd)
    result_dmd, error_dmd = romberg(xrange_dmd, y_data_dmd) 
    # interp_dmd = LinearInterpolation(x_data_dmd, y_data_dmd)
    # result_dmd, error_dmd = quadgk(interp_dmd, minimum(x_data_dmd), maximum(x_data_dmd))

    println("Numerical integration has error: ", max(error_dmd, error_org))
    T1 = 1.0 / abs(result_dmd + result_org) * 0.5
    println("T1 time is: ", T1, ".")
    return T1
end



function process_file(filename, early_truncate::Bool)
    T1 = nothing
    if (filename == "dmd.h5")
        (x_data_dmd, y_data_dmd) = get_kernel_data("dmd.h5", get_ghz_value(filename), early_truncate=early_truncate)
        (x_data_org, y_data_org) = get_kernel_data("orig.h5", get_ghz_value(filename), early_truncate=early_truncate)
        T1 = calculate_T1(x_data_dmd, y_data_dmd, x_data_org, y_data_org)
    else
        (x_data, y_data) = get_kernel_data(filename, get_ghz_value(filename), early_truncate=early_truncate)
        T1 = calculate_T1(x_data, y_data)
    end
    if T1 != nothing
        return T1 
    end
end

function process_directory(base_path, outfile_name, ghz)
    # Get a list of all .h5 files in the directory
    T_dirs = filter(x -> occursin(r"T_*", x), readdir("."))

    T_array = zeros(Float64, length(T_dirs))
    T1_array = zero(T_array)

    for (ii, T_dir) in enumerate(T_dirs)
        file_path = joinpath(base_path, T_dir, "dmd.h5")
        T_array[ii] = get_lmd(T_dir)
        T1_array[ii] = process_file(file_path, ghz)
    end
    # save lmd_array and T1_array in a text file called file_name
    data = hcat(T_array, T1_array)  # Combine arrays horizontally
    header = ["# T"  "T1"]
    writedlm(outfile_name, [header; data], '\t')   # Write as tab-delimited file
end

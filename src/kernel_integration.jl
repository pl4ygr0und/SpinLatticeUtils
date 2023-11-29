using HDF5: h5open
using DelimitedFiles: writedlm, readdlm
# using Interpolations, QuadGK
using Romberg

using AutomaticDocstrings

export get_kernel_data,
    calculate_T1,
    process_file

# util functions
"""
    vec2range(x::Vector)

Transform a evenly spaced vector to a range object. 
The output can be used as an input to the `Romberg` integration module, which restrict its input as a range.
"""
function vec2range(x::Vector)
    xrange = range(minimum(x), maximum(x), length=length(x))
    if collect(xrange) ≈ x
        return xrange
    end
    error("When you want to convert a vector to range, you have input a vector is not evenly spaced!")
end

function get_mask(x_data, y_data, T; early_truncate=false)
    if (T > x_data[end])
        error("Quadratic kernels are oscillatory. Thus you have to provide an origianl data longer than 1 period.")
    else
        dvdval = x_data[end] / T
        fval = floor(Integer, dvdval)
        cval = ceil(Integer, dvdval)
        maskval = isapprox(dvdval, cval; atol=1e-3) ? cval : fval

        return early_truncate ? (x_data .< 1.0 * T) : (x_data .< maskval * T)
    end
end

# main utilties

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

function get_kernel_data(filename)
    if get_extension(filename) == "h5"
        (x_data, y_data) = read_h5_kernel(filename)
    elseif get_extension(filename) == "dat"
        (x_data, y_data) = read_dat_kernel(filename)
    else
        error("You data file isn't 'prop-kernel-eq.dat' or 'dmd.h5'. Check your script! Abort!")
    end

    tol = 1e-6
    devi = abs(y_data[end] / y_data[1])
    if (devi < tol)
        return x_data, y_data
    else
        error("'prop-kernel-eq.dat' not converged. Please use dmd to extend; or you can bump up t_final")
    end
end

function get_kernel_data(filename, ghz; early_truncate=false, mask=true)
    x_data, y_data = get_kernel_data(filename)

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

function compute_k(x_data, y_data)
    xrange = vec2range(x_data)
    result, error = romberg(xrange, y_data)
    return result, error
end

function calculate_T1(
    t_1_gets_0_data,
    K_1_gets_0_data,
    t_0_gets_1_data,
    K_0_gets_1_data
)
    k_1_gets_0, error_1_gets_0 = compute_k(t_1_gets_0_data, K_1_gets_0_data)
    k_0_gets_1, error_0_gets_1 = compute_k(t_0_gets_1_data, K_0_gets_1_data)

    println("="^10, " ", "Numerical integration error: ", "="^10)
    println("Error for k_{1 \\gets 0} is ", error_1_gets_0, ".")
    println("Error for k_{0 \\gets 1} is ", error_0_gets_1, ".")
    k_tot = abs(k_1_gets_0) + abs(k_0_gets_1)
    T1 = 1.0 / k_tot
    println("T1 time is: ", T1, ".")
    println()

    return T1
end

# function calculate_T1(x_data_dmd, y_data_dmd, x_data_org, y_data_org)
#     xrange_org = vec2range(x_data_org)
#     result_org, error_org = romberg(xrange_org, y_data_org) 
#     # interp_org = LinearInterpolation(x_data_org, y_data_org)
#     # result_org, error_org = quadgk(interp_org, minimum(x_data_org), maximum(x_data_org))
# 
#     mask_dmd = (x_data_dmd .>= x_data_org[end])
#     x_data_dmd = x_data_dmd[mask_dmd]
#     y_data_dmd = y_data_dmd[mask_dmd]
#     xrange_dmd = vec2range(x_data_dmd)
#     result_dmd, error_dmd = romberg(xrange_dmd, y_data_dmd) 
#     # interp_dmd = LinearInterpolation(x_data_dmd, y_data_dmd)
#     # result_dmd, error_dmd = quadgk(interp_dmd, minimum(x_data_dmd), maximum(x_data_dmd))
# 
#     println("Numerical integration has error: ", max(error_dmd, error_org))
#     T1 = 1.0 / abs(result_dmd + result_org) * 0.5
#     println("T1 time is: ", T1, ".")
#     return T1
# end


# function process_file(filename, early_truncate::Bool)
#     T1 = nothing
#     if (filename == "dmd.h5")
#         (x_data_dmd, y_data_dmd) = get_kernel_data("dmd.h5", get_ghz_value(filename), early_truncate=early_truncate)
#         (x_data_org, y_data_org) = get_kernel_data("orig.h5", get_ghz_value(filename), early_truncate=early_truncate)
#         T1 = calculate_T1(x_data_dmd, y_data_dmd, x_data_org, y_data_org)
#     else
#         (x_data, y_data) = get_kernel_data(filename, get_ghz_value(filename), early_truncate=early_truncate)
#         T1 = calculate_T1(x_data, y_data)
#     end
#     if T1 != nothing
#         return T1 
#     end
# end
# 
# function process_directory(base_path, outfile_name, ghz)
#     # Get a list of all .h5 files in the directory
#     T_dirs = filter(x -> occursin(r"T_*", x), readdir("."))
# 
#     T_array = zeros(Float64, length(T_dirs))
#     T1_array = zero(T_array)
# 
#     for (ii, T_dir) in enumerate(T_dirs)
#         file_path = joinpath(base_path, T_dir, "dmd.h5")
#         T_array[ii] = get_lmd(T_dir)
#         T1_array[ii] = process_file(file_path, ghz)
#     end
#     # save lmd_array and T1_array in a text file called file_name
#     data = hcat(T_array, T1_array)  # Combine arrays horizontally
#     header = ["# T"  "T1"]
#     writedlm(outfile_name, [header; data], '\t')   # Write as tab-delimited file
# end
# 

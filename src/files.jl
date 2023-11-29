# Conventions about the file systems that is used to organize the data

export get_ghz_value,
    get_T_value,
    get_temperature_value,
    get_λ_value,
    get_lmd_value

"""
    get_extension(filedir)

get the extension of an file
"""
get_extension(filedir) = split(basename(filedir), ".")[end]

"""
    match_int_value(pattern, str)

Use regex expression to match precisely *ONE* integert number. 
Caution: Use this function if you know the matched pattern will yield *only* one integer number.
"""
function match_int_value(pattern, str)
    m = match(pattern, str)
    if m !== nothing
        return int_value = parse(Float64, match(r"\d+", m.match).match)
    else
        error("The string '$(str)' does not match the integer pattern '$(pattern)'.")
    end
end

"""
    match_float_value(pattern, str)

Use regex expression to match precisely *ONE* float number. 
Caution: Use this function if you know the matched pattern will yield *only* one float number.
"""
function match_float_value(pattern, str)
    m = match(pattern, str)
    if m !== nothing
        return float_value = parse(Float64, match(r"\d+.\d+", m.match).match)
    else
        error("The string '$(str)' does not match the float pattern '$(pattern)'.")
    end
end

"""
    get_ghz_value(dir)

from the full directory, get the ghz value
"""
function get_ghz_value(dir)
    abs_dir = abspath(dir) #-- get the abs_dir
    return match_int_value(r"\d+-GHz-equiv", abs_dir)
end

"""
    get_λ_value(dir)

from the full directory, get the λ value
"""
function get_λ_value(dir)
    abs_dir = abspath(dir) #-- get the abs_dir
    return match_float_value(r"lmd_\d+.\d+", abs_dir)
end

"""
    get_temperature_value(dir)

from the full directory, get the temperature value
"""
function get_temperature_value(dir)
    abs_dir = abspath(dir) #-- get the abs_dir
    return match_float_value(r"T_\d+.\d+", abs_dir)
end

"""
    get_T_value(dir)

Alias for get_temperature_value. This exist for old data.
"""
get_T_value(dir) = get_temperature_value(dir)

get_lmd_value(dir) = get_λ_value(dir)

function generate_hash_key(params...)
    # Concatenate all parameters into a tuple
    param_tuple = tuple(params...)
    # Generate hash from the tuple
    return hash(param_tuple)
end

function map_params_to_T1!(mapping, T1, params...)
    hash_key = generate_hash_key(params...)
    param_tuple = tuple(params...)
    mapping[hash_key] = Dict("params" => param_tuple, "T1" => T1)
end

function get_T1_from_params(mapping, params...)
    hash_key = generate_hash_key(params...)
    if haskey(mapping, hash_key)
        if haskey(mapping[hash_key], "T1")
            return mapping[hash_key]["T1"]
        else
            # Raise an error if the "T1" value is not computed for the parameters
            error_msg = "The T1 value for parameter set $params is not computed."
            throw(DomainError(error_msg))
        end
    else
        # Handle the case where the params are not found in the mapping
        error("Parameters not found in the mapping. $(hash_key)")
    end
end

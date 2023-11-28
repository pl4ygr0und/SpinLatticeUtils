module SpinLatticeUtils

using HDF5
using DelimitedFiles
using Romberg
# using Interpolations, QuadGK

using AutomaticDocstrings # add AutomaticDocstrings for a developing code

include("files.jl")
include("SpinLatticeParams.jl")
include("kernel_integration.jl")

end

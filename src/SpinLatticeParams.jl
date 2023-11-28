export SpinLatticeParams

using AutomaticDocstrings

"""
    SpinLatticeParams

A parameter class to hold the input parameters for the spin lattice project.

# Fields:
- `zeeman::Float64`: the Zeeman splitting of the two level
- `coherence::Float64`: the off diagonal coherence of the two level 
- `β::Float64`: the inverse temperature of the bath 
- `λ::Float64`: the system-bath interaction strength
- `Vsph_diag::Float64`: system bath interaction mode: diagonal
- `Vsph_offd::Float64`: system bath interaction mode: off-diagonal 
- `interaction_scheme::String`: the interaction scheme, can only be linear or quad
- `α0::Float64`: quadratic bath mode, the coefficient for the constant  term
- `α1::Float64`: quadratic bath mode, the coefficient for the x_B term
- `α2::Float64`: quadratic bath mode, the coefficient for the x_B^2 term
"""
struct SpinLatticeParams
    zeeman::Float64
    coherence::Float64
    β::Float64
    λ::Float64
    Vsph_diag::Float64
    Vsph_offd::Float64
    interaction_scheme::String
    α0::Float64
    α1::Float64
    α2::Float64
end

"""
    SpinLatticeParams(zeeman::Float64, β::Float64, λ::Float64; coherence::Float64 = 0.0, Vsph_diag::Float64 = 0.0, Vsph_offd::Float64 = 1.0, interaction_scheme::String = linear, α0::Float64 = 0.0, α1::Float64 = 0.0, α2::Float64 = 0.0)

The constructor for the SpinLatticeParams struct

# Arguments:
- `zeeman`: see SpinLatticeParams
- `β`: see SpinLatticeParams
- `λ`: see SpinLatticeParams
- `coherence`: defaults to 0. No intrinsic coherence
- `Vsph_diag`: defaults to 0; i.e., system-bath interaction don't reorganalize the energy level
- `Vsph_offd`: defaults to 1.0; meaning system-bath coupling will cause unit off diagonal transition
- `interaction_scheme`: defaults to 'linear'
- `α0`: defaults to 0.0, since assuming linear coupling.
- `α1`: defaults to 0.0, since assuming linear coupling.
- `α2`: defaults to 0.0, since assuming linear coupling.
"""
function SpinLatticeParams(
    zeeman::Float64,
    β::Float64,
    λ::Float64;
    coherence::Float64=0.0,
    Vsph_diag::Float64=0.0,
    Vsph_offd::Float64=1.0,
    interaction_scheme::String="linear",
    α0::Float64=0.0,
    α1::Float64=0.0,
    α2::Float64=0.0)
    if interaction_scheme != "linear" && interaction_scheme != "quadratic"
        throw(ArgumentError("Please input interaction_scheme `linear` or `quadratic`."))
    end

    if interaction_scheme == "linear"
        return SpinLatticeParams(zeeman, coherence, β, λ, Vsph_diag, Vsph_offd, interaction_scheme, 0.0, 0.0, 0.0)
    else
        return SpinLatticeParams(zeeman, coherence, β, λ, Vsph_diag, Vsph_offd, interaction_scheme, α0, α1, α2)
    end
end

"""
    disp_mat(mat; tab = true, tabn = 1)

Helper function to display a matrix.

# Arguments:
- `mat`: the matrix to be displayed.
- `tab`: add tab before the the matrix or not.
- `tabn`: the number of tabs to be added. Defaults to 1.
"""
disp_mat(mat; tab=true, tabn=1) = join(["\t"^tabn * s for s in split(string(mat), '\n')], "\n")

"""
    get_H(slp::SpinLatticeParams)

Get the system Hamiltonian from the slp struct.
"""
get_H(slp::SpinLatticeParams) = [slp.zeeman slp.coherence; slp.coherence -slp.zeeman]

"""
    get_Q(slp::SpinLatticeParams)

Get the general system mode to interact with the system
"""
get_Q(slp::SpinLatticeParams) = [0.0 slp.Vsph_offd; slp.Vsph_offd slp.Vsph_diag]

"""
    linear_str(slp::SpinLatticeParams)

Sumarize the most important information about a linear SpinLatticeParams in to a string. To be printed out to the console.
"""
function linear_str(slp::SpinLatticeParams)
    if slp.interaction_scheme != "linear"
        error("This is the method generates a string for a quadratic SpinLatticeParams. But your `interaction_scheme` is $(slp.interaction_scheme).")
    end
    s = "\n"
    s *= "Spin-Lattice interaction model with linear system-bath coupling.\n"
    s *= "- The system Hamiltonian:\n"
    s *= "$(disp_mat(get_H(slp)))\n"
    s *= "- The heat bath temperature:\n"
    s *= "\t$(1/slp.β) (inverse temperature is $(slp.β))\n"
    s *= "- The system mode interact with a boson bath:\n"
    s *= "$(disp_mat(get_Q(slp)))\n"
    s *= "- The interaction strength:\n"
    s *= "\t$(slp.λ)"

    return s
end

"""
    quad_str(slp::SpinLatticeParams)

Sumarize the most important information about a quadratic SpinLatticeParams in to a string. To be printed out to the console.
"""
function quad_str(slp::SpinLatticeParams)
    if slp.interaction_scheme != "quadratic"
        error("This is the method generates a string for a quadratic SpinLatticeParams. But your `interaction_scheme` is $(slp.interaction_scheme).")
    end
    s = "\n"
    s *= "Spin-Lattice interaction model with quadratic system-bath coupling.\n"
    s *= "- The system Hamiltonian:\n"
    s *= "$(disp_mat(get_H(slp)))\n"
    s *= "- The heat bath temperature:\n"
    s *= "\t$(1/slp.β) (inverse temperature is $(slp.β))\n"
    s *= "- The system mode interact with a boson bath:\n"
    s *= "$(disp_mat(get_Q(slp)))\n"
    s *= "- The bosonic bath interaction operators are\n"
    s *= "\tF = $(slp.α0) + $(slp.α1)*x_B + $(slp.α2)*x_B^2\n"
    s *= "- The interaction strength:\n"
    s *= "\t$(slp.λ)"

    return s
end

function Base.show(io::IO, slp::SpinLatticeParams)
    if slp.interaction_scheme == "linear"
        println(io, linear_str(slp))
    elseif slp.interaction_scheme == "quadratic"
        println(io, quad_str(slp))
    end
end

# small helper functions about the sipn lattice system
dimless_energy_2_ghz(E::Number) = 100.0 * E

ghz_2_dimless_energy(ghz::Number) = ghz / 100.0

"""
    ghz_2_period(ghz::Number)

Calculate the peroid from the Ghz number. Note for the 100.0, this is I use 
    1 dimensionless energy unit = 100 GHz 
to make some comparison with the experiment people.
"""
ghz_2_period(ghz::Number) = π * 2.0 * 100. / ghz

"""
    ghz_2_period(ghz::AbstractString)

Calculate the peroid from the GHZ string.
"""
ghz_2_period(ghz::AbstractString) = get_period(parse(Float64, ghz))
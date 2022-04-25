

"""
    AbstractWavefunction

Abstract supertype for all wave functions.
"""
abstract type AbstractWavefunction end


"""
    RBM <: AbstractWavefunction

Represent a Restricted Boltzmann Machine wave function.
"""
struct RBM{T<:Number} <: AbstractWavefunction
    a::Vector{T}        # Bias vector of visible units
    b::Vector{T}        # Bias vector of hidden units
    W::Matrix{T}        # Weights connecting visible and hidden units
    θ::Vector{T}        # Look-up table of effective angles
end


function RBM{T}(Nσ::Int, Nh::Int) where {T<:Number}
    @assert(Nσ > 0 && Nh > 0)

    return RBM(ones(T,Nσ), ones(T,Nh), ones(T,Nh,Nσ), ones(T,Nh))
end


"""Return the number of visible units `Nσ`"""
get_Nσ(rbm::RBM) = length(rbm.a)

"""Return the number of hidden units `Nh`"""
get_Nh(rbm::RBM) = length(rbm.b)


"""
    log_psi(rbm::RBM, σ::Vector{Int})

Return the log of the RBM wave function amplitude for state σ.
"""
function log_psi(rbm::RBM, σ::Vector{Int})

    logΨ = dot(rbm.a, σ)

    for i = 1:get_Nh(rbm)
        θ = rbm.b[i] + dot(rbm.W[i,:], σ)
        logΨ += log(cosh(θ))
    end

    return logΨ
end

"""
    log_psi_over_psi(rbm::RBM, σ::Vector{Int}, flips::Vector{Int})

Computes the logarithm of Ψ(σ′)/Ψ(σ) where σ′ is state σ with sites flipped
according to the entries in `flips`.
"""
function log_psi_over_psi(rbm::RBM, σ::Vector{Int}, flips::Vector{Int})

    if length(flips) == 0
        return Complex(0.0,0.0)
    end

    logΨoverΨ = 0.0 + 0.0im

    # changes due to bias
    for j in flips
        logΨoverΨ -= 2.0 * rbm.a[j] * σ[j]
    end

    for i = 1:get_Nh(rbm)
        # get effective angle from look-up table
        θ = rbm.θ[i]

        # changes due weights
        for j in flips
            θ -= 2.0 * rbm.W[i,j] * σ[j]
        end

        logΨoverΨ += log(cosh(θ)) - log(cosh(rbm.θ[i]))
    end

    return logΨoverΨ

end


"""
    psi_over_psi(rbm::RBM, σ::Vector{Int}, flips::Vector{Int})

Computes Ψ(σ′)/Ψ(σ) where σ′ is state σ with sites flipped
according to the entries in `flips`.
"""
psi_over_psi(rbm::RBM, σ::Vector{Int}, flips::Vector{Int}) = exp(log_psi_over_psi(rbm, σ, flips))


function update_θ(rbm::RBM, σ::Vector{Int}, flips::Vector{Int})

    if length(flips) == 0
        return
    end

    for i = 1:get_Nh(rbm)
        for j in flips
            rbm.θ[i] -= 2.0 * rbm.W[i,j] * σ[j]
        end
    end

    return

end

module UKMUQ

using Reexport

@reexport using Distributions

using QuasiMonteCarlo

# export structs/methods

export evaluate
export probability_of_failure
export sample
export Model

export AbstractMonteCarlo
export MonteCarlo
export SobolSampling

abstract type AbstractMonteCarlo end

struct MonteCarlo <: AbstractMonteCarlo
    n::Int
end

struct SobolSampling <: AbstractMonteCarlo
    n::Int
end

function sample(rvs::Vector{T} where {T<:UnivariateDistribution}, sim::MonteCarlo)
    m = length(rvs)
    A = zeros(sim.n, m)

    for i in 1:m
        A[:, i] = rand(rvs[i], sim.n)
    end

    return A
end

function sample(rvs::Vector{T} where {T<:UnivariateDistribution}, sim::SobolSampling)


    m = length(rvs)
    x = QuasiMonteCarlo.sample(sim.n, zeros(m), ones(m), SobolSample())'

    for i in 1:m
        x[:, i] = quantile.(rvs[i], x[:, i])
    end

    return Matrix(x)

end

struct Model
    f::Function
end

function evaluate(m::Model, x::Matrix)
    return m.f(x)
end

function probability_of_failure(rvs::Vector{T} where {T<:UnivariateDistribution}, m::Model, performance::Function, sim::AbstractMonteCarlo)
    x = sample(rvs, sim)
    y = evaluate(m, x)

    g = performance(y)

    pf = sum(g .< 0) / sim.n

    cov = (sqrt(pf - pf^2) / sim.n) / pf

    return pf, cov
end

end # module UKMUQ

using UKMUQ

function displacement(x::Matrix)
    I = (0.12 .* x[:, 1] .^ 3) ./ 12

    w = (x[:, 4] .* 9.81 .* 0.12 .* x[:, 1] .* 1.8^4) ./ (8 .* x[:, 2] .* I)
    w = w .+ (x[:, 3] .* 1.8^3) ./ (3 .* x[:, 2] .* I)
    return w
end

m = Model(displacement)

rvs = [Normal(0.24, 0.01), LogNormal(23.0, 0.16), LogNormal(8.5, 0.08), LogNormal(6.37, 0.23)]

n = 1000000
pf, cov = probability_of_failure(rvs, m, y -> 0.01 .- y, n)

println("Probability of failure: $pf ($cov)")

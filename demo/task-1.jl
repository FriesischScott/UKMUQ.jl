using UKMUQ
using Plots

rvs = [Normal(), Normal(1.0, 1.0), Uniform(-2, 2)]
n = 10000

samples = sample(rvs, n)

histogram(samples[:, 1], alpha=0.5, bins=100)
histogram!(samples[:, 2], alpha=0.5, bins=100)
histogram!(samples[:, 3], alpha=0.5, bins=100)

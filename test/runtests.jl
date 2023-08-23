import SD4SOLPS
using Plots
using Test

function make_test_profile(x; sym=0.9657, hwid=0.05421, offset=0.0953, pedestal=3.58901, alpha=0.005)
    z = (sym .- x) ./ hwid
    amplitude = (pedestal - offset) / 2.0
    b = offset + amplitude
    # modified tanh profile formula from:
    # Groebner, et al., Nucl. Fusion 41, 1789 (2001) 10.1088/0029-5515/41/12/306
    return b .+ amplitude .* ((1 .+ alpha .* z) .* exp.(z) .- exp.(.-z)) ./ (exp.(z) .+ exp.(.-z))
end

function plot_test_profiles()
    edge_rho = Array(LinRange(0.88, 1.0, 18))
    edge_quantity = make_test_profile(edge_rho)
    output_rho = Array(LinRange(0, 1.0, 201))
    output_quantity = SD4SOLPS.extrapolate_core(edge_rho, edge_quantity, output_rho)
    Plots.plot(output_rho, output_quantity)
    Plots.plot!(edge_rho, edge_quantity, marker='.')
end

@testset "core_profile_extension" begin
    edge_rho = Array(LinRange(0.88, 1.0, 18))
    edge_quantity = make_test_profile(edge_rho)
    output_rho = Array(LinRange(0, 1.0, 201))
    output_quantity = SD4SOLPS.extrapolate_core(edge_rho, edge_quantity, output_rho)
    @test length(output_quantity) == length(output_rho)
end


import SD4SOLPS
import SOLPS2IMAS
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

# function plot_test_profiles()
#     edge_rho = Array(LinRange(0.88, 1.0, 18))
#     edge_quantity = make_test_profile(edge_rho)
#     output_rho = Array(LinRange(0, 1.0, 201))
#     output_quantity = SD4SOLPS.extrapolate_core(edge_rho, edge_quantity, output_rho)
#     Plots.plot(output_rho, output_quantity)
#     Plots.plot!(edge_rho, edge_quantity, marker='.')
# end

@testset "core_profile_extension" begin
    edge_rho = Array(LinRange(0.88, 1.0, 18))
    edge_quantity = make_test_profile(edge_rho)
    output_rho = Array(LinRange(0, 1.0, 201))
    output_quantity = SD4SOLPS.extrapolate_core(edge_rho, edge_quantity, output_rho)
    @test length(output_quantity) == length(output_rho)
end

@testset "utilities" begin
    sample_path = splitdir(pathof(SOLPS2IMAS))[1] * "/../samples/"
    file_list = SD4SOLPS.find_files_in_allowed_folders(sample_path, eqdsk_file="thereisntoneyet")
    @test length(file_list) == 4
    b2fgmtry, b2time, gridspec, eqdsk = file_list
    @test length(b2fgmtry) > 10
    @test endswith(b2fgmtry, "b2fgmtry")
    @test length(b2time) > 10
    @test endswith(b2time, "b2time.nc")
    @test length(gridspec) > 10
    @test endswith(gridspec, "gridspacedesc.yml")
end

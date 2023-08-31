import SD4SOLPS
import SOLPS2IMAS
import OMAS
import EFIT
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
    file_list = SD4SOLPS.find_files_in_allowed_folders(
        sample_path, eqdsk_file="thereisntoneyet"
    )
    @test length(file_list) == 4
    b2fgmtry, b2time, gridspec, eqdsk = file_list
    @test length(b2fgmtry) > 10
    @test endswith(b2fgmtry, "b2fgmtry")
    @test length(b2time) > 10
    @test endswith(b2time, "b2time.nc")
    @test length(gridspec) > 10
    @test endswith(gridspec, "gridspacedesc.yml")
end

@testset "geqdsk_to_imas" begin
    sample_files = (splitdir(pathof(SD4SOLPS))[1] * "/../sample/") .* [
        "g184833.03600", "geqdsk_iter_small_sample"
    ]
    tslice = 1
    for sample_file in sample_files
        dd = OMAS.dd()
        SD4SOLPS.geqdsk_to_imas(sample_file, dd, time_index=tslice)
        eqt = dd.equilibrium.time_slice[tslice]

        # global
        gq = eqt.global_quantities
        @test gq.magnetic_axis.r > 0
        @test dd.equilibrium.vacuum_toroidal_field.r0 > 0

        # 1d
        p1 = eqt.profiles_1d
        nprof = length(p1.psi)
        @test nprof > 10
        @test p1.psi[1] == gq.psi_axis
        @test p1.psi[end] == gq.psi_boundary
        @test length(p1.f) == nprof
        @test length(p1.pressure) == nprof
        @test length(p1.rho_tor_norm) == nprof

        # 2d
        p2 = eqt.profiles_2d[1]
        @test length(p2.grid.dim1) > 10
        @test length(p2.grid.dim2) > 10
        @test size(p2.psi) == (length(p2.grid.dim1), length(p2.grid.dim2))

        # derived
        @test gq.q_axis == p1.q[1]

        # boundary
        @test length(eqt.boundary.outline.r) == length(eqt.boundary.outline.z)

        # wall
        limiter = dd.wall.description_2d[1].limiter
        @test length(limiter.unit[1].outline.r) > 10
        @test length(limiter.unit[1].outline.r) == length(limiter.unit[1].outline.z)
    end
end

@testset "preparation" begin
    eqdsk_file = "geqdsk_iter_small_sample"
    sample_paths = [
        splitdir(pathof(SOLPS2IMAS))[1] * "/../samples/",
        splitdir(pathof(SD4SOLPS))[1] * "/../sample/",
    ]
    dd = SD4SOLPS.preparation(eqdsk_file, sample_paths...)
    p2 = dd.equilibrium.time_slice[1].profiles_2d[1]
    psirz = p2.psi
    r = p2.grid.dim1
    z = p2.grid.dim2
    @test size(psirz) == (length(r), length(z))
end

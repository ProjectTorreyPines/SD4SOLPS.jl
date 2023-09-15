import SD4SOLPS
import SOLPS2IMAS
import OMAS
import EFIT
using Plots
using Test

"""
    make_test_profile()

Makes a modified tanh profile based on the formula from:
Groebner, et al., Nucl. Fusion 41, 1789 (2001) 10.1088/0029-5515/41/12/306

x: position basis. The output will be a function of this coordinate.
    Could be psi_N or R or Z or whatever, as long as sym and hwid are
    specified in terms of the same coordinate.
    Defaults are given in terms of psi_N.
sym: location of the tanh symmetry point
hwid: half width of the tanh
offset: value at the foot of the tanh / SOL value
pedestal: value at the top of the tanh / pedestal value
alpha: interior slope
"""
function make_test_profile(x; sym=0.9657, hwid=0.05421, offset=0.0953, pedestal=3.58901, alpha=0.005)
    z = (sym .- x) ./ hwid
    amplitude = (pedestal - offset) / 2.0
    b = offset + amplitude
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
    # Just the basic profile extrapolator
    edge_rho = Array(LinRange(0.88, 1.0, 18))
    edge_quantity = make_test_profile(edge_rho)
    output_rho = Array(LinRange(0, 1.0, 201))
    output_quantity = SD4SOLPS.extrapolate_core(edge_rho, edge_quantity, output_rho)
    @test length(output_quantity) == length(output_rho)

    # The full workflow
    sample_path = splitdir(pathof(SD4SOLPS))[1] * "/../sample/ITER_Lore_2296_00000/extended_output"
    sample_path2 = splitdir(pathof(SD4SOLPS))[1] * "/../sample/ITER_Lore_2296_00000/baserun"
    sample_path3 = splitdir(pathof(SD4SOLPS))[1] * "/../sample/ITER_Lore_2296_00000/run_restart"
    sample_path4 = splitdir(pathof(SOLPS2IMAS))[1] * "/../samples"

    # Requires dvc pull before the full samples will be loaded
    # Sorry, the minimal samples are too minimal for this one.
    file_list = SD4SOLPS.find_files_in_allowed_folders(
        sample_path, sample_path2, sample_path3, sample_path4,
        eqdsk_file="thereisntoneyet",
        allow_reduced_versions=false,
    )
    b2fgmtry, b2time, b2mn, gridspec, eqdsk = file_list
    eqdsk = splitdir(pathof(SD4SOLPS))[1] * "/../sample/ITER_Lore_2296_00000/EQDSK/Baseline2008-li0.70.x4.mod3.eqdsk"
    println(b2fgmtry)
    println(b2time)
    println(b2mn)
    println(gridspec)
    println(eqdsk)
    # If these files don't exist, complete the DVC sample setup and try again
    @test isfile(b2fgmtry)
    @test isfile(b2time)
    @test isfile(b2mn)
    @test isfile(gridspec)
    @test isfile(eqdsk)
    dd = SOLPS2IMAS.solps2imas(b2fgmtry, b2time, gridspec, b2mn)    
    SD4SOLPS.geqdsk_to_imas(eqdsk, dd)
    SD4SOLPS.fill_in_extrapolated_core_profile(dd, "electrons.density")
end

@testset "utilities" begin
    # Test for finding files in allowed folders
    sample_path = splitdir(pathof(SOLPS2IMAS))[1] * "/../samples/"
    file_list = SD4SOLPS.find_files_in_allowed_folders(
        sample_path, eqdsk_file="thereisntoneyet", allow_reduced_versions
    )
    @test length(file_list) == 5
    b2fgmtry, b2time, b2mn, gridspec, eqdsk = file_list
    @test length(b2fgmtry) > 10
    @test endswith(b2fgmtry, "b2fgmtry_red")
    @test length(b2time) > 10
    @test endswith(b2time, "b2time_red.nc")
    @test length(b2mn) > 10
    @test endswith(b2mn, "b2mn.dat")
    @test length(gridspec) > 10
    @test endswith(gridspec, "gridspacedesc.yml")

    # Test for sweeping 1D core profiles into 2D R,Z (or anyway evaluating them at any R,Z location)
    dd = OMAS.dd()
    eqdsk_file = splitdir(pathof(SD4SOLPS))[1] * "/../sample/geqdsk_iter_small_sample"
    SD4SOLPS.geqdsk_to_imas(eqdsk_file, dd)
    quantity = "electrons.density"
    prof_time_idx = eq_time_idx = 1
    resize!(dd.core_profiles.profiles_1d, prof_time_idx)
    n = 101
    rho_n = Array(LinRange(0, 1.0, n))
    resize!(dd.core_profiles.profiles_1d[prof_time_idx].grid.rho_tor_norm, n)
    resize!(dd.core_profiles.profiles_1d[prof_time_idx].electrons.density, n)
    dd.core_profiles.profiles_1d[prof_time_idx].grid.rho_tor_norm = rho_n
    dd.core_profiles.profiles_1d[prof_time_idx].electrons.density = make_test_profile(rho_n)
    r_mag_axis = dd.equilibrium.time_slice[1].global_quantities.magnetic_axis.r
    z_mag_axis = dd.equilibrium.time_slice[1].global_quantities.magnetic_axis.z
    rg = r_mag_axis .* [0.85, 0.90, 0.95, 1.0, 1.05, 1.10, 1.15]
    zg = z_mag_axis .+ (r_mag_axis .* [-0.05, -0.025, 0, 0.025, 0.05])
    points = collect(Base.Iterators.product(rg, zg))
    r = getfield.(points, 1)
    z = getfield.(points, 2)
    density_on_grid = SD4SOLPS.core_profile_2d(dd, prof_time_idx, eq_time_idx, quantity, r, z)
    @test size(density_on_grid) == (length(rg), length(zg))
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

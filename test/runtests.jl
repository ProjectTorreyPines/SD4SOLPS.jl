using SD4SOLPS: SD4SOLPS
using SOLPS2IMAS: SOLPS2IMAS
using OMAS: OMAS
using EFIT: EFIT
using Plots
using Test
using Unitful: Unitful
using Interpolations: Interpolations
using ArgParse: ArgParse

function parse_commandline()
    s = ArgParse.ArgParseSettings(; description="Run tests. Default is all tests.")

    ArgParse.add_arg_table!(s,
        ["--lightweight_utilities"],
        Dict(:help => "Test only lightweight utilities",
            :action => :store_true),
        ["--actuator"],
        Dict(:help => "Test only actuator model",
            :action => :store_true),
        ["--core_profile_extension"],
        Dict(:help => "Test only core profile extension",
            :action => :store_true),
        ["--edge_profile_extension"],
        Dict(:help => "Test only edge profile extension",
            :action => :store_true),
        ["--heavy_utilities"],
        Dict(:help => "Test only heavy utilities",
            :action => :store_true),
        ["--repair_eq"],
        Dict(:help => "Test only repair_eq",
            :action => :store_true),
        ["--geqdsk_to_imas"],
        Dict(:help => "Test only geqdsk_to_imas",
            :action => :store_true),
        ["--preparation"],
        Dict(:help => "Test only preparation",
            :action => :store_true),
    )
    args = ArgParse.parse_args(s)
    if !any(values(args)) # If no flags are set, run all tests
        for k ∈ keys(args)
            args[k] = true
        end
    end
    return args
end
args = parse_commandline()

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
function make_test_profile(
    x;
    sym=0.9657,
    hwid=0.05421,
    offset=0.0953,
    pedestal=3.58901,
    alpha=0.005,
)
    z = (sym .- x) ./ hwid
    amplitude = (pedestal - offset) / 2.0
    b = offset + amplitude
    return b .+
           amplitude .* ((1 .+ alpha .* z) .* exp.(z) .- exp.(.-z)) ./
           (exp.(z) .+ exp.(.-z))
end

# function plot_test_profiles()
#     edge_rho = Array(LinRange(0.88, 1.0, 18))
#     edge_quantity = make_test_profile(edge_rho)
#     output_rho = Array(LinRange(0, 1.0, 201))
#     output_quantity = SD4SOLPS.extrapolate_core(edge_rho, edge_quantity, output_rho)
#     Plots.plot(output_rho, output_quantity)
#     Plots.plot!(edge_rho, edge_quantity, marker='.')
# end

function define_default_sample_set()
    sample_path =
        splitdir(pathof(SD4SOLPS))[1] *
        "/../sample/ITER_Lore_2296_00000/extended_output"
    sample_path2 =
        splitdir(pathof(SD4SOLPS))[1] * "/../sample/ITER_Lore_2296_00000/baserun"
    sample_path3 =
        splitdir(pathof(SD4SOLPS))[1] * "/../sample/ITER_Lore_2296_00000/run_restart"
    sample_path4 = splitdir(pathof(SOLPS2IMAS))[1] * "/../samples"

    # Requires dvc pull before the full samples will be loaded
    # Sorry, the minimal samples are too minimal for this one.
    file_list = SD4SOLPS.find_files_in_allowed_folders(
        sample_path, sample_path2, sample_path3, sample_path4;
        eqdsk_file="thereisntoneyet",
        allow_reduced_versions=false,
    )
    b2fgmtry, b2time, b2mn, gridspec, eqdsk = file_list
    eqdsk =
        splitdir(pathof(SD4SOLPS))[1] *
        "/../sample/ITER_Lore_2296_00000/EQDSK/Baseline2008-li0.70.x4.mod2.eqdsk"
    return b2fgmtry, b2time, b2mn, gridspec, eqdsk
end

if args["lightweight_utilities"]
    @testset "lightweight_utilities" begin
        # Gas unit converter
        flow_tls = 40.63 * Unitful.Torr * Unitful.L / Unitful.s
        flow_pam3 = SD4SOLPS.gas_unit_converter(flow_tls, "torr L s^-1", "Pa m^3 s^-1")
        flow_pam3_no_unitful =
            SD4SOLPS.gas_unit_converter(flow_tls.val, "torr L s^-1", "Pa m^3 s^-1")
        @test flow_pam3.val > 0.0
        @test flow_pam3_no_unitful ==
              Unitful.uconvert(Unitful.Pa * Unitful.m^3 / Unitful.s, flow_pam3).val

        flow_molecules1 = SD4SOLPS.gas_unit_converter(
            flow_tls,
            "torr L s^-1",
            "molecules s^-1";
            temperature=293.15 * Unitful.K,
        )
        flow_molecules2 = SD4SOLPS.gas_unit_converter(
            flow_tls,
            "torr L s^-1",
            "molecules s^-1";
            temperature=300.0 * Unitful.K,
        )
        flow_molecules3 = SD4SOLPS.gas_unit_converter(
            flow_tls.val,
            "torr L s^-1",
            "molecules s^-1";
            temperature=300.0,
        )
        @test flow_molecules1 > flow_molecules2
        @test (Unitful.upreferred(flow_molecules2)).val == flow_molecules3
    end
end

if args["actuator"]
    @testset "actuator" begin
        t = collect(LinRange(0, 2.0, 1001))
        cmd = (t .> 1.0) * 1.55 + (t .> 1.5) * 0.93 + (t .> 1.8) * 0.25

        instant_flow_function = SD4SOLPS.model_gas_valve("instant")
        instant_flow = instant_flow_function(t, cmd)
        @test length(instant_flow) == length(cmd)

        simple_flow_function = SD4SOLPS.model_gas_valve("simple")
        simple_flow = simple_flow_function(t, cmd)
        @test length(simple_flow) == length(cmd)
    end
end

if args["core_profile_extension"]
    @testset "core_profile_extension" begin
        # Just the basic profile extrapolator ------------------
        edge_rho = Array(LinRange(0.88, 1.0, 18))
        edge_quantity = make_test_profile(edge_rho)
        output_rho = Array(LinRange(0, 1.0, 201))
        output_quantity = SD4SOLPS.extrapolate_core(edge_rho, edge_quantity, output_rho)
        @test length(output_quantity) == length(output_rho)

        # The full workflow --------------------------------------
        # Setup sample DD
        b2fgmtry, b2time, b2mn, gridspec, eqdsk = define_default_sample_set()
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
        rho = dd.equilibrium.time_slice[1].profiles_1d.rho_tor_norm

        if !SD4SOLPS.check_rho_1d(dd; time_slice=1)
            SD4SOLPS.add_rho_to_equilibrium!(dd)
            rho = dd.equilibrium.time_slice[1].profiles_1d.rho_tor_norm
            println("Repaired missing rho for core profile test")
        end
        # Sample is ready

        # Test settings
        quantity_name = "electrons.density"
        test_slice_idx = 1

        # Do it
        SD4SOLPS.fill_in_extrapolated_core_profile!(dd, quantity_name)

        # Inspect results
        @test length(dd.core_profiles.profiles_1d) > 0
        core_prof = dd.core_profiles.profiles_1d[test_slice_idx]
        tags = split(quantity_name, ".")
        quantity = dd.core_profiles.profiles_1d[test_slice_idx]
        for tag ∈ tags
            quantity = getproperty(quantity, Symbol(tag))
        end
        rho_core = dd.core_profiles.profiles_1d[test_slice_idx].grid.rho_tor_norm
        @test length(quantity) > 0
        @test length(quantity) == length(rho_core)
    end
end

if args["edge_profile_extension"]
    @testset "edge_profile_extension" begin
        # Test for getting mesh spacing
        b2fgmtry, b2time, b2mn, gridspec, eqdsk = define_default_sample_set()
        dd = SOLPS2IMAS.solps2imas(b2fgmtry, b2time, gridspec, b2mn)
        SD4SOLPS.geqdsk_to_imas(eqdsk, dd)
        dpsin = SD4SOLPS.mesh_psi_spacing(dd)
        @test dpsin > 0.0

        n_edge = 47
        n_outer_prof = 13
        quantity_edge =
            [collect(LinRange(2.0, 75.0, 7)); collect(LinRange(75.0, 98.3, 17))]
        quantity_edge = [quantity_edge; reverse(quantity_edge)[2:end]]
        gradient_edge =
            [collect(LinRange(-1, -100, 7)); collect(LinRange(-100, -350, 17))]
        gradient_edge = [gradient_edge; reverse(gradient_edge)[2:end]]
        psin_out = collect(LinRange(1.0, 1.25, n_outer_prof + 1))[2:end]
    end
end

if args["heavy_utilities"]
    @testset "heavy_utilities" begin
        # Test for finding files in allowed folders
        sample_path = splitdir(pathof(SOLPS2IMAS))[1] * "/../samples"
        file_list = SD4SOLPS.find_files_in_allowed_folders(
            sample_path; eqdsk_file="thereisntoneyet", allow_reduced_versions=true,
        )
        @test length(file_list) == 5
        b2fgmtry, b2time, b2mn, gridspec, eqdsk = file_list
        @test length(b2fgmtry) > 10
        @test endswith(b2fgmtry, "b2fgmtry_red") | endswith(b2fgmtry, "b2fgmtry")
        @test length(b2time) > 10
        @test endswith(b2time, "b2time_red.nc") | endswith(b2time, "b2time.nc")
        @test length(b2mn) > 10
        @test endswith(b2mn, "b2mn.dat")
        @test length(gridspec) > 10
        @test endswith(gridspec, "gridspacedesc.yml")

        # Test for sweeping 1D core profiles into 2D R,Z
        # (or anyway evaluating them at any R,Z location)
        dd = OMAS.dd()
        eqdsk_file =
            splitdir(pathof(SD4SOLPS))[1] * "/../sample/geqdsk_iter_small_sample"
        SD4SOLPS.geqdsk_to_imas(eqdsk_file, dd)
        quantity = "electrons.density"
        prof_time_idx = eq_time_idx = 1
        resize!(dd.core_profiles.profiles_1d, prof_time_idx)
        n = 101
        rho_n = Array(LinRange(0, 1.0, n))
        resize!(dd.core_profiles.profiles_1d[prof_time_idx].grid.rho_tor_norm, n)
        resize!(dd.core_profiles.profiles_1d[prof_time_idx].electrons.density, n)
        dd.core_profiles.profiles_1d[prof_time_idx].grid.rho_tor_norm = rho_n
        dd.core_profiles.profiles_1d[prof_time_idx].electrons.density =
            make_test_profile(rho_n)
        r_mag_axis = dd.equilibrium.time_slice[1].global_quantities.magnetic_axis.r
        z_mag_axis = dd.equilibrium.time_slice[1].global_quantities.magnetic_axis.z
        rg = r_mag_axis .* [0.85, 0.90, 0.95, 1.0, 1.05, 1.10, 1.15]
        zg = z_mag_axis .+ (r_mag_axis .* [-0.05, -0.025, 0, 0.025, 0.05])
        points = collect(Base.Iterators.product(rg, zg))
        r = getfield.(points, 1)
        z = getfield.(points, 2)
        if !SD4SOLPS.check_rho_1d(dd; time_slice=eq_time_idx)
            SD4SOLPS.add_rho_to_equilibrium!(dd)
            println("DD was repaired (rho added) for core 2d utility test")
        end
        density_on_grid =
            SD4SOLPS.core_profile_2d(dd, prof_time_idx, eq_time_idx, quantity, r, z)
        @test size(density_on_grid) == (length(rg), length(zg))
    end
end

if args["repair_eq"]
    @testset "repair_eq" begin
        # Prepare sample
        dd = OMAS.dd()
        eqdsk = splitdir(pathof(SD4SOLPS))[1] * "/../sample/geqdsk_iter_small_sample"
        SD4SOLPS.geqdsk_to_imas(eqdsk, dd)
        # Make sure rho is missing
        nt = length(dd.equilibrium.time_slice)
        for it ∈ 1:nt
            eqt = dd.equilibrium.time_slice[it]
            eqt.profiles_1d.rho_tor_norm = Vector{Float64}()
        end
        # Add rho
        SD4SOLPS.add_rho_to_equilibrium!(dd)
        # Check
        rho = dd.equilibrium.time_slice[1].profiles_1d.rho_tor_norm
        println(rho)
        @test length(rho) > 10
        @test maximum(rho) > 0
    end
end

if args["geqdsk_to_imas"]
    @testset "geqdsk_to_imas" begin
        sample_files =
            (splitdir(pathof(SD4SOLPS))[1] * "/../sample/") .* [
                "g184833.03600", "geqdsk_iter_small_sample",
            ]
        append!(
            sample_files,
            ["$(@__DIR__)/../sample/ITER_Lore_2296_00000/EQDSK/g002296.00200"],
        )
        tslice = 1
        for sample_file ∈ sample_files
            println(sample_file)
            dd = OMAS.dd()
            SD4SOLPS.geqdsk_to_imas(sample_file, dd; time_index=tslice)
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

            # boundary
            r_bry = eqt.boundary.outline.r
            z_bry = eqt.boundary.outline.z
            @test length(r_bry) == length(z_bry)

            # 2d
            # Did the R and Z (dim1 and dim2) coordinates get written properly?
            p2 = eqt.profiles_2d[1]
            r_eq = p2.grid.dim1
            z_eq = p2.grid.dim2
            @test length(r_eq) > 10
            @test minimum(r_eq) > 0  # R should be dim1
            @test minimum(z_eq) < 0
            @test length(z_eq) > 10
            # Does psi match R and Z?
            println(size(p2.psi), (length(r_eq), length(z_eq)))
            @test size(p2.psi) == (length(r_eq), length(z_eq))
            # Does the psi grid look like it's transposed the right way?
            # Many equilibrium reconstructions have equal numbers of cells in R and Z,
            # so transposing the psi map incorrectly would not be detected by checking
            # its array dimensions. It's also possible to incorrectly associate the psi
            # map with R and Z (dim1 and dim2). So what recourse does that leave us?
            # Well, the points on the boundary outline should be at psi_N = 1.0 to
            # within some small tolerance.
            psin2d = (p2.psi .- gq.psi_axis) ./ (gq.psi_boundary - gq.psi_axis)
            tolerance = 2.0e-3  # It's not always a high res contour so cut some slack
            psin_bry =
                Interpolations.LinearInterpolation((r_eq, z_eq), psin2d).(r_bry, z_bry)
            @test maximum(psin_bry) < (1.0 + tolerance)
            @test minimum(psin_bry) > (1.0 - tolerance)

            # derived
            @test gq.q_axis == p1.q[1]

            # wall
            limiter = dd.wall.description_2d[1].limiter
            @test length(limiter.unit[1].outline.r) > 10
            @test length(limiter.unit[1].outline.r) == length(limiter.unit[1].outline.z)
        end
    end
end

if args["preparation"]
    @testset "preparation" begin
        eqdsk_file = "geqdsk_iter_small_sample"
        sample_paths = [
            splitdir(pathof(SD4SOLPS))[1] * "/../sample",
            splitdir(pathof(SOLPS2IMAS))[1] * "/../samples",
        ]
        core_method = "simple"
        edge_method = "simple"
        filename = splitdir(pathof(SD4SOLPS))[1] * "/../sd_input_data"
        output_format = "json"
        dd = SD4SOLPS.preparation(
            eqdsk_file,
            sample_paths...;
            core_method=core_method,
            edge_method=edge_method,
            filename=filename,
            output_format=output_format,
        )
        out_file = filename * "." * output_format
        p2 = dd.equilibrium.time_slice[1].profiles_2d[1]
        psirz = p2.psi
        r = p2.grid.dim1
        z = p2.grid.dim2
        @test size(psirz) == (length(r), length(z))
        println(out_file)
        @test isfile(out_file)
    end
end

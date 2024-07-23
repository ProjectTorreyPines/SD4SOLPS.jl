module SD4SOLPS

using IMASDD: IMASDD
using SOLPS2IMAS: SOLPS2IMAS
using EFIT: EFIT
using Interpolations: Interpolations

export find_files_in_allowed_folders, geqdsk_to_imas!, preparation

include("$(@__DIR__)/supersize_profile.jl")
include("$(@__DIR__)/repair_eq.jl")
include("$(@__DIR__)/unit_utils.jl")

"""
    find_files_in_allowed_folders(
        input_dirs::String...;
        eqdsk_file::String,
        recursive::Bool=true,
    )

Searches a list of allowed folders for a set of filenames that will provide information
about the SOLPS case. Returns a list of filenames with complete paths.

Example:

```julia
SD4SOLPS.find_files_in_allowed_folders(
    "<your samples folder>/D3D_Ma_184833_03600";
    eqdsk_file="g184833.03600",
)
```
"""
function find_files_in_allowed_folders(
    input_dirs::String...;
    eqdsk_file::String,
    recursive::Bool=true,
)::Vector{String}
    files = ["b2fgmtry", "b2time.nc", "b2mn.dat", eqdsk_file]
    reduced_files =
        ["b2fgmtry_red", "b2time_red.nc", "b2mn.dat", eqdsk_file]
    output_files = fill("", length(files))
    if recursive
        dirs = []
        for dir ∈ input_dirs
            dirs =
                append!(dirs, [subdir[1] for subdir ∈ [item for item ∈ walkdir(dir)]])
        end
    else
        dirs = input_dirs
    end
    for i ∈ eachindex(files)
        for dir ∈ dirs
            file = dir * "/" * files[i]
            reduced_file = dir * "/" * reduced_files[i]
            if isfile(file)
                output_files[i] = file
                break
            elseif isfile(reduced_file)
                output_files[i] = reduced_file
                break
            end
        end
    end
    return output_files
end

"""
    geqdsk_to_imas!(
        eqdsk_file::String,
        dd::IMASDD.dd;
        set_time::Union{Nothing, Float64}=nothing,
        time_index::Int=1,
    )

Transfers the equilibrium reconstruction from an EFIT-style gEQDSK file into
the IMAS DD structure.
"""
function geqdsk_to_imas!(
    eqdsk_file::String,
    dd::IMASDD.dd;
    set_time::Union{Nothing, Float64}=nothing,
    time_index::Int=1,
)
    # https://github.com/JuliaFusion/EFIT.jl/blob/master/src/io.jl
    g = EFIT.readg(eqdsk_file; set_time=set_time)
    gfilename = split(eqdsk_file, "/")[end]
    # Copying ideas from OMFIT: omfit/omfit_classes/omfit_eqdsk.py / to_omas()
    eq = dd.equilibrium
    if IMASDD.ismissing(eq, :time)
        eq.time = Array{Float64}(undef, time_index)
    end
    eq.time[time_index] = g.time
    if length(eq.time_slice) < time_index
        resize!(eq.time_slice, time_index)
    end
    eqt = eq.time_slice[time_index]
    eqt.time = g.time

    source_for_summary = "gEQDSK file $gfilename loaded during SD4SOLPS workflow."

    # 0D
    gq = eqt.global_quantities
    gq.magnetic_axis.r = g.rmaxis
    gq.magnetic_axis.z = g.zmaxis
    gq.psi_axis = g.simag
    gq.psi_boundary = g.sibry
    gq.ip = g.current
    eq.vacuum_toroidal_field.r0 = g.rcentr
    b0 = Array{Float64}(undef, time_index)
    b0[time_index] = g.bcentr
    eq.vacuum_toroidal_field.b0 = b0

    if IMASDD.ismissing(dd.summary, :time)
        dd.summary.time = Array{Float64}(undef, time_index)
    end
    dd.summary.time[time_index] = g.time
    ip = Array{Float64}(undef, time_index)
    ip[time_index] = g.current
    dd.summary.global_quantities.ip.value = ip
    dd.summary.global_quantities.r0.value = g.rcentr
    dd.summary.global_quantities.b0.value = b0
    summarize = ["ip", "r0", "b0"]

    # 1D
    p1 = eqt.profiles_1d
    nprof = length(g.pres)
    psi = collect(LinRange(g.simag, g.sibry, nprof))
    p1.psi = psi
    p1.f = g.fpol
    p1.pressure = g.pres
    p1.f_df_dpsi = g.ffprim
    p1.dpressure_dpsi = g.pprime
    p1.q = g.qpsi
    if hasproperty(g, :rhovn)
        # rhovn is not in the original EFIT.jl but is added on a branch
        p1.rho_tor_norm = g.rhovn
    end

    # 2D
    resize!(eqt.profiles_2d, 1)
    p2 = eqt.profiles_2d[1]
    p2.grid.dim1 = collect(g.r)
    p2.grid.dim2 = collect(g.z)
    p2.psi = g.psirz  # Not sure if transpose is correct
    # missing j_tor = pcurrt

    # Derived
    psin1d = (psi .- g.simag) ./ (g.sibry - g.simag)
    gq.magnetic_axis.b_field_tor = g.bcentr * g.rcentr / g.rmaxis
    gq.q_axis = g.qpsi[1]
    gq.q_95 = Interpolations.linear_interpolation(psin1d, g.qpsi)(0.95)
    qmin_idx = argmin(abs.(g.qpsi))
    gq.q_min.value = g.qpsi[qmin_idx]
    if hasproperty(g, :rhovn)
        gq.q_min.rho_tor_norm = g.rhovn[qmin_idx]
    end

    dd.summary.global_quantities.q_95.value = Array{Float64}(undef, time_index)
    dd.summary.global_quantities.q_95.value[time_index] = gq.q_95
    summarize = [summarize; "q_95"]

    # X-points
    xrs, xzs, xpsins, xseps = EFIT.x_points(g; within_limiter_only=false)
    if length(xrs) > 0
        bx = eqt.boundary.x_point
        resize!(bx, length(xrs))
        for i ∈ eachindex(xrs)
            bx[i].r = xrs[i]
            bx[i].z = xzs[i]
        end
        nprim = sum(xseps .== 1)
        if nprim > 0
            bsx = eqt.boundary_separatrix.x_point
            resize!(bsx, nprim)
            xrprim = xrs[xseps.==1]
            xzprim = xzs[xseps.==1]
            for i ∈ nprim
                bsx[i].r = xrprim[i]
                bsx[i].z = xzprim[i]
            end
        end
        nsec = sum(xseps .== 2)
        if nsec > 0
            bssx = eqt.boundary_secondary_separatrix.x_point
            resize!(bssx, nsec)
            xrsec = xrs[xseps.==2]
            xzsec = xzs[xseps.==2]
            for i ∈ nsec
                bssx[i].r = xrsec[i]
                bssx[i].z = xzsec[i]
            end
        end
    end

    # Boundary / LCFS
    eqt.boundary.outline.r = g.rbbbs
    eqt.boundary.outline.z = g.zbbbs

    # Wall
    resize!(dd.wall.description_2d, 1)
    limiter = dd.wall.description_2d[1].limiter
    limiter.type.name = "first wall"
    limiter.type.index = 0
    limiter.type.description = "first wall"
    resize!(limiter.unit, 1)
    limiter.unit[1].outline.r = g.rlim
    limiter.unit[1].outline.z = g.zlim

    for s ∈ summarize
        getproperty(dd.summary.global_quantities, Symbol(s)).source = source_for_summary
    end

    return
end

"""
    preparation(
        eqdsk_file::String,
        dirs::String...;
        core_method::String="simple",
        filename::String="sd_input_data",
        output_format::String="json",
        eqdsk_set_time::Union{Nothing, Float64}=nothing,
        eq_time_index::Int64=1,
    )::IMASDD.dd

Gathers SOLPS and EFIT files and loads them into IMAS structure. Extrapolates
profiles as needed to get a complete picture.
"""
function preparation(
    eqdsk_file::String,
    dirs::Vector{String};
    core_method::String="simple",
    filename::String="sd_input_data",
    output_format::String="json",
    eqdsk_set_time::Union{Nothing, Float64}=nothing,
    eq_time_index::Int64=1,
)::IMASDD.dd
    b2fgmtry, b2time, b2mn, eqdsk =
        find_files_in_allowed_folders(dirs...; eqdsk_file=eqdsk_file)
    println("Found source files:")
    println("    b2fgmtry = ", b2fgmtry)
    println("    b2time = ", b2time)
    println("    b2mn.dat = ", b2mn)
    println("    eqdsk = ", eqdsk)

    dd = SOLPS2IMAS.solps2imas(b2fgmtry, b2time; b2mn=b2mn)
    geqdsk_to_imas!(eqdsk, dd; set_time=eqdsk_set_time, time_index=eq_time_index)
    # Repairs
    add_rho_to_equilibrium!(dd)  # Doesn't do anything if rho is valid
    println("Loaded input data into IMAS DD")

    core_profiles = ["electrons.density", "electrons.temperature"]
    extrapolated_core_profiles = []
    for core_profile ∈ core_profiles
        tags = split(core_profile, ".")
        parent = dd.edge_profiles.ggd[1]
        for tag ∈ tags[1:end-1]
            parent = getproperty(parent, Symbol(tag))
        end
        qty = getproperty(parent, Symbol(tags[end]), core_profile)
        if length(qty) > 0
            fill_in_extrapolated_core_profile!(dd, core_profile; method=core_method)
            append!(extrapolated_core_profiles, [core_profile])
            println("    Extrapolated $core_profile into the core.")
        else
            println(
                "  > Warning: quantity $core_profile was not usable and was not extrapolated into the core.",
            )
        end
    end
    # ... more profiles here as they become available in b2time
    println(
        "Extrapolated $(length(extrapolated_core_profiles))/$(length(core_profiles)) core profiles.",
    )

    cached_mesh_extension!(dd, eqdsk_file, b2fgmtry)
    fill_in_extrapolated_edge_profile!(dd, "electrons.density"; method=core_method)
    # ... more profiles here
    println("Extrapolated edge profiles (but not really (placeholder only))")

    if output_format == "json"
        IMASDD.imas2json(dd, filename * ".json")
    else
        throw(ArgumentError(string("Unrecognized output format: ", output_format)))
    end
    return dd
end

end # module SD4SOLPS

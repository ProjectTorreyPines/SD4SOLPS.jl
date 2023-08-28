module SD4SOLPS

import OMAS
import SOLPS2IMAS
import EFIT
import Interpolations
#import GGDUtils

export find_files_in_allowed_folders
export geqdsk_to_imas

include("$(@__DIR__)/supersize_profile.jl")

greet() = print("Hello World!")

"""
    find_files_in_allowed_folders()

Searches a list of allowed folders for a set of filenames that will provide information
about the SOLPS case. Returns a list of filenames with complete paths.

Example:
SD4SOLPS.find_files_in_allowed_folders("<your samples folder>/D3D_Ma_184833_03600", eqdsk_file="g184833.03600")
"""
function find_files_in_allowed_folders(input_dirs...; eqdsk_file, recursive=true)
    files = ["b2fgmtry", "b2time.nc", "gridspacedesc.yml", eqdsk_file]
    println(files)
    output_files = fill("", length(files))
    if recursive
        dirs = []
        for dir in input_dirs
            dirs = append!(dirs, [subdir[1] for subdir in [item for item in walkdir(dir)]])
        end
    else
        dirs = input_dirs
    end
    for i in 1:length(files)
        for dir in dirs
            file = dir * "/" * files[i]
            if isfile(file)
                output_files[i] = file
                break
            end
        end
    end
    return output_files
end

function geqdsk_to_imas(eqdsk_file, dd; time_index=1)
    """
        geqdsk_to_imas()

    Transfers the equilibrium reconstruction in an EFIT-style gEQDSK file into
    the IMAS DD structure.
    """
    # https://github.com/JuliaFusion/EFIT.jl/blob/master/src/io.jl
    g = EFIT.readg(eqdsk_file)
    # Copying ideas from OMFIT: omfit/omfit_classes/omfit_eqdsk.py / to_omas()
    eq = dd.equilibrium
    resize!(eq.time_slice, 1)
    eqt = eq.time_slice[time_index]

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
    if hasproperty(g, :rhovn)  # rhovn is not in the original EFIT.jl but is added on a branch
        p1.rho_tor_norm = g.rhovn
    end

    # 2D
    resize!(eqt.profiles_2d, 1)
    p2 = eqt.profiles_2d[1]
    p2.grid.dim1 = collect(g.r)
    p2.grid.dim2 = collect(g.z)
    p2.psi = Matrix(transpose(g.psirz))  # Not sure if transpose is correct
    # missing j_tor = pcurrt

    # Derived
    psin1d = (psi .- g.simag) ./ (g.sibry - g.simag)
    gq.magnetic_axis.b_field_tor = g.bcentr * g.rcentr / g.rmaxis
    gq.q_axis = g.qpsi[1]
    gq.q_95 = Interpolations.LinearInterpolation(psin1d, g.qpsi)(0.95)
    qmin_idx = argmin(abs.(g.qpsi))
    gq.q_min.value = g.qpsi[qmin_idx]
    if hasproperty(g, :rhovn)
        gq.q_min.rho_tor_norm = g.rhovn[qmin_idx]
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

end

"""
    preparation()

Gathers SOLPS and EFIT files and loads them into IMAS structure. Extrapolates
profiles as needed to get a complete picture.
"""
function preparation(eqdsk_file, dirs...)
    b2fgmtry, b2time, gridspec, eqdsk = find_files_in_allowed_folders(dirs, eqdsk_file=eqdsk_file)
    dd = SOLPS2IMAS.solps2imas(b2gmtry, b2output, gsdesc)
    geqdsk_to_imas(eqdsk_file, dd)
end

end # module SD4SOLPS

module SD4SOLPS

import OMAS
import SOLPS2IMAS
#import GGDUtils

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

"""
    preparation()

Gathers SOLPS and EFIT files and loads them into IMAS structure. Extrapolates
profiles as needed to get a complete picture.
"""
function preparation(eqdsk_file, dirs...)
    b2fgmtry, b2time, gridspec, eqdsk = find_files_in_allowed_folders(dirs, eqdsk_file=eqdsk_file)
    dd = SOLPS2IMAS.solps2imas(b2gmtry, b2output, gsdesc)

    # https://github.com/JuliaFusion/EFIT.jl/blob/master/src/io.jl
    gfile = EFIT.readg(eqdsk)

end

end # module SD4SOLPS

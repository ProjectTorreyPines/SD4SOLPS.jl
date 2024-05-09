var documenterSearchIndex = {"docs":
[{"location":"#SD4SOLPS.jl","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"Pages = [\"index.md\"]\nDepth = 3","category":"page"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"This repository serves as the top most workflow manager with helpful utilities to use other repositories in this project.","category":"page"},{"location":"#Documentation-of-other-repositories-in-this-project","page":"SD4SOLPS.jl","title":"Documentation of other repositories in this project","text":"","category":"section"},{"location":"#[GGDUtils.jl](https://projecttorreypines.github.io/GGDUtils.jl/stable)","page":"SD4SOLPS.jl","title":"GGDUtils.jl","text":"","category":"section"},{"location":"#[SOLPS2IMAS.jl](https://projecttorreypines.github.io/SOLPS2IMAS.jl/stable)","page":"SD4SOLPS.jl","title":"SOLPS2IMAS.jl","text":"","category":"section"},{"location":"#[SynthDiag.jl](https://projecttorreypines.github.io/SynthDiag.jl/stable)","page":"SD4SOLPS.jl","title":"SynthDiag.jl","text":"","category":"section"},{"location":"#Installation","page":"SD4SOLPS.jl","title":"Installation","text":"","category":"section"},{"location":"#Using-make:","page":"SD4SOLPS.jl","title":"Using make:","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"After cloning this repo, check the make menu:","category":"page"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"SD4SOLPS.jl % make help\nHelp Menu\n\nmake env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories\nmake env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones\nmake clean: Deletes Project.toml and Manifest.toml for a fresh start","category":"page"},{"location":"#make-r","page":"SD4SOLPS.jl","title":"make r","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"This option creates local copies of required private repositories at the same level as current repository and uses them in develop mode to create a Manifest.toml","category":"page"},{"location":"#make-u","page":"SD4SOLPS.jl","title":"make u","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"This option uses url of required private repositories to create a static Manifest.toml attached to current master branches of these repositories.","category":"page"},{"location":"#make-clean","page":"SD4SOLPS.jl","title":"make clean","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"Deletes Manifest.toml so that environment can be recreated, to update or change the last used method.","category":"page"},{"location":"#Using-Julia-REPL-and-installing-using-Github-url","page":"SD4SOLPS.jl","title":"Using Julia REPL and installing using Github url","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"julia> using Pkg;\njulia> Pkg.add(; url=\"https://github.com/ProjectTorreyPines/IMASDD.jl.git\");\njulia> Pkg.add(; url=\"https://github.com/ProjectTorreyPines/GGDUtils.jl.git\");\njulia> Pkg.add(; url=\"https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl.git\");\njulia> Pkg.add(; url=\"https://github.com/JuliaFusion/EFIT.jl.git\");\njulia> Pkg.add(; url=\"https://github.com/ProjectTorreyPines/SD4SOLPS.jl.git\");\njulia> Pkg.instantiate()","category":"page"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"You might have to use ssh url instead of https. In that case, replace https://github.com with git@github.com:.","category":"page"},{"location":"#Top-file-handling-functions","page":"SD4SOLPS.jl","title":"Top file handling functions","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"find_files_in_allowed_folders\ngeqdsk_to_imas!\npreparation","category":"page"},{"location":"#SD4SOLPS.find_files_in_allowed_folders","page":"SD4SOLPS.jl","title":"SD4SOLPS.find_files_in_allowed_folders","text":"find_files_in_allowed_folders(\n    input_dirs::String...;\n    eqdsk_file::String,\n    recursive::Bool=true,\n)\n\nSearches a list of allowed folders for a set of filenames that will provide information about the SOLPS case. Returns a list of filenames with complete paths.\n\nExample:\n\nSD4SOLPS.find_files_in_allowed_folders(\n    \"<your samples folder>/D3D_Ma_184833_03600\";\n    eqdsk_file=\"g184833.03600\",\n)\n\n\n\n\n\n","category":"function"},{"location":"#SD4SOLPS.geqdsk_to_imas!","page":"SD4SOLPS.jl","title":"SD4SOLPS.geqdsk_to_imas!","text":"geqdsk_to_imas!(\n    eqdsk_file::String,\n    dd::IMASDD.dd;\n    set_time::Union{Nothing, Float64}=nothing,\n    time_index::Int=1,\n)\n\nTransfers the equilibrium reconstruction from an EFIT-style gEQDSK file into the IMAS DD structure.\n\n\n\n\n\n","category":"function"},{"location":"#SD4SOLPS.preparation","page":"SD4SOLPS.jl","title":"SD4SOLPS.preparation","text":"preparation(\n    eqdsk_file::String,\n    dirs::String...;\n    core_method::String=\"simple\",\n    filename::String=\"sd_input_data\",\n    output_format::String=\"json\",\n    eqdsk_set_time::Union{Nothing, Float64}=nothing,\n    eq_time_index::Int64=1,\n)::IMASDD.dd\n\nGathers SOLPS and EFIT files and loads them into IMAS structure. Extrapolates profiles as needed to get a complete picture.\n\n\n\n\n\n","category":"function"},{"location":"#Repairing/filling-out-partial-equilibrium-files","page":"SD4SOLPS.jl","title":"Repairing/filling out partial equilibrium files","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"Tools for repairing/filling out partial equilibrium files.","category":"page"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"Some of the added fields may not be totally accurate, so it is recommended to use this tool mainly for test cases, as a utility. For a real equilibrium, problems should be fixed properly.","category":"page"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"add_rho_to_equilibrium!\ncheck_rho_1d","category":"page"},{"location":"#SD4SOLPS.add_rho_to_equilibrium!","page":"SD4SOLPS.jl","title":"SD4SOLPS.add_rho_to_equilibrium!","text":"function add_rho_to_equilibrium(dd:IMASDD.dd)\n\nAdds equilibrium rho profile to the DD\n\n\n\n\n\n","category":"function"},{"location":"#SD4SOLPS.check_rho_1d","page":"SD4SOLPS.jl","title":"SD4SOLPS.check_rho_1d","text":"check_rho_1d(\n    dd::IMASDD.dd;\n    time_slice::Int64=1,\n    throw_on_fail::Bool=false,\n)::Bool\n\nChecks to see if rho exists and is valid in the equilibrium 1d profiles\n\n\n\n\n\n","category":"function"},{"location":"#Extrapolations","page":"SD4SOLPS.jl","title":"Extrapolations","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"Utilities for extrapolating profiles","category":"page"},{"location":"#Core-profile-extrapolations","page":"SD4SOLPS.jl","title":"Core profile extrapolations","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"extrapolate_core\nfill_in_extrapolated_core_profile!","category":"page"},{"location":"#SD4SOLPS.extrapolate_core","page":"SD4SOLPS.jl","title":"SD4SOLPS.extrapolate_core","text":"extrapolate_core(\n    edge_rho::Vector{Float64},\n    edge_quantity::Vector{Float64},\n    rho_output::Vector{Float64},\n)::Vector{Float64}\n\nFunction for assuming a core profile when given edge profile data.\n\nConcept:\n\nvalue and derivative should be continuous when joining real data\nderivative at magnetic axis is known to be 0 when making profile vs. rho, by the definition of rho\nderivative probably does something fancier between the pedestal and axis than just linear interpolation, so add an extra point in there\nthere's a joint between the steep pedestal and the shallow core that needs an extra knot to manage it properly\nafter making up a very simple gradient profile out of a few line segments, integrate it to get the profile of the quantity in question\n\n\n\n\n\n","category":"function"},{"location":"#SD4SOLPS.fill_in_extrapolated_core_profile!","page":"SD4SOLPS.jl","title":"SD4SOLPS.fill_in_extrapolated_core_profile!","text":"fill_in_extrapolated_core_profile!(\n    dd::IMASDD.dd,\n    quantity_name::String;\n    method::String=\"simple\",\n    eq_time_idx::Int64=1,\n    eq_profiles_2d_idx::Int64=1,\n    grid_ggd_idx::Int64=1,\n    space_idx::Int64=1,\n)\n\nThis function accepts a DD that should be populated with equilibrium and edge_profiles as well as a request for a quantity to extrapolate into the core. It then maps edge_profiles data to rho, calls the function that performs the extrapolation (which is not a simple linear extrapolation but has some trickery to attempt to make a somewhat convincing profile shape), and writes the result to core_profiles. This involves a bunch of interpolations and stuff.\n\nInput arguments:\n\ndd: an IMAS data dictionary\nquantity_name: the name of a quantity in edge_profiles.profiles_2d and core_profiles.profiles_1d, such as \"electrons.density\"\nmethod: Extrapolation method.\neq_time_idx: index of the equilibrium time slice to use. For a typical SOLPS run, the SOLPS mesh will be based on the equilibrium reconstruction at a single time, so the DD associated with the SOLPS run only needs one equilibrium time slice to be loaded. However, one could combine the complete equilibrium time series with the SOLPS run and then have to specify which slice of the equilibrium corresponds to the SOLPS mesh.\neq_profiles_2d_idx: index of the profiles_2D in equilibrium time_slice.\ngrid_ggd_idx: index of the grid_ggd to use. For a typical SOLPS run, the SOLPS grid is fixed, so this index defaults to 1. But in future, if a time varying grid is used, then this index will need to be specified.\nspace_idx: index of the space to use. For a typical SOLPS run, there will be only one space so this index will mostly remain at 1.\ncell_subset_idx: index of the subset of cells to use for the extrapolation. The default is 5, which is the subset of all cells. If edge_profiles data is instead present for a different subset, for instance, -5, which are b2.5 cells only, then this index should be set to -5.\n\n\n\n\n\n","category":"function"},{"location":"#Edge-profiles-extrapolations","page":"SD4SOLPS.jl","title":"Edge profiles extrapolations","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"These functions have not been fully tested and/or supported yet.","category":"page"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"mesh_psi_spacing\ncached_mesh_extension!","category":"page"},{"location":"#SD4SOLPS.mesh_psi_spacing","page":"SD4SOLPS.jl","title":"SD4SOLPS.mesh_psi_spacing","text":"mesh_psi_spacing(\n    dd::IMASDD.dd;\n    eq_time_idx::Int64=1,\n    eq_profiles_2d_idx::Int64=1,\n    grid_ggd_idx::Int64=1,\n    space_idx::Int64=1,\n    avoid_guard_cell::Bool=true,\n    spacing_rule=\"mean\",\n)\n\nInspects the mesh to see how far apart faces are in psi_N. Requires that GGD and equilibrium are populated.\n\nInput Arguments:\n\ndd: a data dictionary instance with required data loaded into it\neq_time_idx: index of the equilibrium time slice to use. For a typical SOLPS run, the SOLPS mesh will be based on the equilibrium reconstruction at a single time, so the DD associated with the SOLPS run only needs one equilibrium time slice to be loaded. However, one could combine the complete equilibrium time series with the SOLPS run and then have to specify which slice of the equilibrium corresponds to the SOLPS mesh.\neq_profiles_2d_idx: index of the profiles_2D in equilibrium time_slice.\ngrid_ggd_idx: index of the grid_ggd to use. For a typical SOLPS run, the SOLPS grid is fixed, so this index defaults to 1. But in future, if a time varying grid is used, then this index will need to be specified.\nspace_idx: index of the space to use. For a typical SOLPS run, there will be only one space so this index will mostly remain at 1.\navoid_guard_cell: assume that the last cell is a guard cell so take end-2 and end-1 instead of end and end-1\nspacing_rule: \"edge\" or \"mean\" to make spacing of new cells (in psi_N) be the same as the spacing at the edge of the mesh, or the same as the average spacing\n\n\n\n\n\n","category":"function"},{"location":"#SD4SOLPS.cached_mesh_extension!","page":"SD4SOLPS.jl","title":"SD4SOLPS.cached_mesh_extension!","text":"cached_mesh_extension!(\n    dd::IMASDD.dd,\n    eqdsk_file::String,\n    b2fgmtry::String;\n    eq_time_idx::Int64=1,\n    eq_profiles_2d_idx::Int64=1,\n    grid_ggd_idx::Int64=1,\n    space_idx::Int64=1,\n    clear_cache::Bool=false,\n)::String\n\nAdds an extended mesh to a data dictionary, possibly from a cached result.\n\nInput Arguments:\n\ndd: The data dictionary. It will be modified in place.\neqdsk_file: the name of the EQDSK file that was used to get equilibrium data in the dd.\nb2fgmtry: the name of the SOLPS geometry file that was used to get GGD info in edge_profiles in the dd.\neq_time_idx: Index of the time slice in equilibrium\neq_profiles_2d_idx: Index of the 2D profile set in equilibrium (there is usually only one)\ngrid_ggd_idx: Index of the grid_ggd set in edge_profiles\nspace_idx: Index of the space\nclear_cache: delete any existing cache file (for use in testing)\n\n\n\n\n\n","category":"function"},{"location":"#Unit-conversion-utilities","page":"SD4SOLPS.jl","title":"Unit conversion utilities","text":"","category":"section"},{"location":"","page":"SD4SOLPS.jl","title":"SD4SOLPS.jl","text":"gas_unit_converter","category":"page"},{"location":"#SD4SOLPS.gas_unit_converter","page":"SD4SOLPS.jl","title":"SD4SOLPS.gas_unit_converter","text":"gas_unit_converter(\n    value_in::Float64,\n    units_in::String,\n    units_out::String;\n    species::String=\"H\",\n    temperature::Float64=293.15,\n)\n\nConverts gas flows between different units. Uses ideal gas law to convert between Pressure * volume type flows / quantities and count / current types of units. There is a version that accepts floats in and outputs floats, and another that deals in Unitful quantities.\n\n\n\n\n\ngas_unit_converter(\n    value_in::Unitful.Quantity,\n    units_in::String,\n    units_out::String;\n    species::String=\"H\",\n    temperature=293.15 * Unitful.K,\n)\n\nConverts gas flows between different units. Uses ideal gas law to convert between Pressure * volume type flows / quantities and count / current types of units. This is the Unitful version.\n\nOutput will be unitful, but the units are not simplified automatically. You can perform operations such as\n\n(output |> Unitful.upreferred).val\nUnitful.uconvert(Unitful.whatever, output).val\n\nto handle simplification or conversion of units.\n\nAlthough this function pretends torr L s^-1 and Pa m^3 s^-1 are different, use of Unitful should cause them to behave the same way as long as you simplify or convert units at the end. This means that you can use other pressure*volume type gas units and call them torr L s^-1 and the script will deal with them up to having messy units in the output.\n\n\n\n\n\n","category":"function"}]
}

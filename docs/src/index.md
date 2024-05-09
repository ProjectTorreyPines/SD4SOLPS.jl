
# SD4SOLPS.jl 

```@contents
Pages = ["index.md"]
Depth = 3
```

This repository serves as the top most workflow manager with helpful utilities to use other repositories in this project.

## Documentation of other repositories in this project

### [GGDUtils.jl](https://projecttorreypines.github.io/GGDUtils.jl/stable)

### [SOLPS2IMAS.jl](https://projecttorreypines.github.io/SOLPS2IMAS.jl/stable)

### [SynthDiag.jl](https://projecttorreypines.github.io/SynthDiag.jl/stable)

## Installation

### Using make:
After cloning this repo, check the make menu:
```
SD4SOLPS.jl % make help
Help Menu

make env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories
make env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones
make clean: Deletes Project.toml and Manifest.toml for a fresh start
```

#### make r
This option creates local copies of required private repositories at the same level as current repository and uses them in develop mode to create a Manifest.toml

#### make u
This option uses url of required private repositories to create a static Manifest.toml attached to current master branches of these repositories.

#### make clean
Deletes Manifest.toml so that environment can be recreated, to update or change the last used method.

### Using Julia REPL and installing using Github url

Or, in julia REPL:
```julia
julia> using Pkg;
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/IMASDD.jl.git");
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/GGDUtils.jl.git");
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl.git");
julia> Pkg.add(; url="https://github.com/JuliaFusion/EFIT.jl.git");
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/SD4SOLPS.jl.git");
julia> Pkg.instantiate()
```

## Top file handling functions

```@docs
find_files_in_allowed_folders
geqdsk_to_imas!
preparation
```

## Repairing/filling out partial equilibrium files

Tools for repairing/filling out partial equilibrium files.

Some of the added fields may not be totally accurate, so it is recommended to
use this tool mainly for test cases, as a utility. For a real equilibrium,
problems should be fixed properly.

```@docs
add_rho_to_equilibrium!
check_rho_1d
```

## Extrapolations

Utilities for extrapolating profiles

### Core profile extrapolations

```@docs
extrapolate_core
fill_in_extrapolated_core_profile!
```

### Edge profiles extrapolations

These functions have not been fully tested and/or supported yet.

```@docs
mesh_psi_spacing
cached_mesh_extension!
```

## Unit conversion utilities

```@docs
gas_unit_converter
```
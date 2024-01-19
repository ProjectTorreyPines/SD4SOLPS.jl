#!/bin/bash

echo "Generating Project.toml and Manifest.toml"
OMAS_url="git@github.com:ProjectTorreyPines/OMAS.jl.git"
EFIT_url="git@github.com:JuliaFusion/EFIT.jl.git"
julia --project=. -e 'using Pkg; Pkg.add(url="'$OMAS_url'", rev="master"); Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.add(url="'$EFIT_url'", rev="master"); Pkg.instantiate()'
for repo in GGDUtils SOLPS2IMAS SynthDiag SD4SOLPS; do
    repo_url="git@github.com:ProjectTorreyPines/"$repo".jl.git"
    julia --project=. -e 'using Pkg; Pkg.add(url="'$repo_url'", rev="master"); Pkg.instantiate()'
done
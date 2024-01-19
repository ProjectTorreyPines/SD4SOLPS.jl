#!/bin/bash

is_repo=$(git rev-parse --is-inside-work-tree 2>/dev/null)
if [[ ${is_repo} == "true" ]]
then
    SD4SOLPS_path="$(git rev-parse --show-toplevel)"
    project_path="$(dirname ${SD4SOLPS_path})"
else
    SD4SOLPS_path=""
    project_path="$(dirname $(pwd))"
fi
echo "Cloning git repositories to "$project_path
orig_path=$(pwd)
cd $project_path
if [[ -z $SD4SOLPS_path ]]
then
    git clone git@github.com:ProjectTorreyPines/SD4SOLPS.jl.git
    SD4SOLPS_path=$project_path"/SD4SOLPS.jl"
fi
# in a for loop clone SynthDiag.jl SOLPS2IMAS.jl OMAS.jl GGDUtils.jl
for repo in SynthDiag.jl SOLPS2IMAS.jl OMAS.jl GGDUtils.jl
do
    if [[ ! -d $repo ]]
    then
        git clone "git@github.com:ProjectTorreyPines/"$repo".git"
    fi
done
cd $orig_path

echo "Generating Project.toml and Manifest.toml"
OMAS_path=$project_path"/OMAS.jl"
EFIT_url="git@github.com:JuliaFusion/EFIT.jl.git"
julia --project=. -e 'using Pkg; Pkg.develop(path="'$OMAS_path'"); Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.add(url="'$EFIT_url'", rev="master"); Pkg.instantiate()'
for repo in GGDUtils SOLPS2IMAS SynthDiag SD4SOLPS; do
    repo_path=$project_path"/"$repo".jl"
    julia --project=. -e 'using Pkg; Pkg.develop(path="'${repo_path}'"); Pkg.instantiate()'
done
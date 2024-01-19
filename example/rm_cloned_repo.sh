#!/bin/bash

echo "Deleting any cloned repo"
is_repo=$(git rev-parse --is-inside-work-tree 2>/dev/null)
if [[ ${is_repo} == "true" ]]
then
    SD4SOLPS_path="$(git rev-parse --show-toplevel)"
    project_path="$(dirname ${SD4SOLPS_path})"
else
    SD4SOLPS_path=""
    project_path="$(dirname $(pwd))"
fi

echo "Removing extra cloned git repositories from "$project_path

if [[ -z $SD4SOLPS_path ]]
then
    # Check if $project_path"/SD4SOLPS.jl" exists
    if [[ -d $project_path"/SD4SOLPS.jl" ]]
    then
        echo "Removing "$project_path"/SD4SOLPS.jl"
        rm -rf $project_path"/SD4SOLPS.jl"
    fi
fi
for repo in SynthDiag.jl SOLPS2IMAS.jl OMAS.jl GGDUtils.jl
do
    if [[ -d $project_path"/"$repo ]]
    then
        echo "Removing "$project_path"/"$repo
        rm -rf $project_path"/"$repo
    fi
done
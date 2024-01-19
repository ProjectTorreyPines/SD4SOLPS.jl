# SD4SOLPS
Synthetic diagnostic workflow manager for use with SOLPS models

This repository is the top level layer for managing a workflow for calculating
the outputs of synthetic diagnostics attached to a SOLPS model.

Steps:
1) Load SOLPS into IMAS DD format if not already
2) Load equilibrium (that the SOLPS mesh was based on) into IMAS DD format
3) Make assumptions to extend profiles into the core and far SOL, if needed
4) Run synthetic diagnostic models and record output


## Installation

### Using make file

#### Option 1: Download the [example directory](https://github.com/ProjectTorreyPines/SD4SOLPS.jl/tree/master/example) in this repo:

```bash
% make help 
Help Menu

make env_with_cloned_repo: Creates a Julia environment with the cloned repositories
make env_with_git_url: Creates a Julia environment with the git urls without creating local clones
make clean: Deletes Project.toml and Manifest.toml for a fresh start
make Clean: Deletes Project.toml, Manifest.toml, and any cloned repositories for a fresh start
```
##### Environment with cloned repos

`make env_with_cloned_repo` will clone all required git repos one level up from your current working directory and use them to define the environment files that will be generated in your working directory. This way, in future you can do `git pull` or change branches to update your working environment.

##### Static environment with git url

`make env_with_git_url` will not create local git clones for you, but it will still clone master branches of all the required repos in julia package management directory and will create a static environment for you that will remain tied to the version of packages when this environment was created.

##### Clean up

`make clean` or `make Clean` can be used clean up environment files and redo the environment generation.

#### Option 2: Clone this repo

```bash
% git clone git@github.com:ProjectTorreyPines/SD4SOLPS.jl.git
% cd example
$ make help
Help Menu

make env_with_cloned_repo: Creates a Julia environment with the cloned repositories
make env_with_git_url: Creates a Julia environment with the git urls without creating local clones
make clean: Deletes Project.toml and Manifest.toml for a fresh start
make Clean: Deletes Project.toml, Manifest.toml, and any cloned repositories for a fresh start
```

Further options are same as above except for the difference that in case of cloning local copies of repos, they will be kept on same level as where you cloned SD4SOLPS.jl repo.

### Manual Install

Refer to the instructions on this [wiki page](https://github.com/ProjectTorreyPines/SD4SOLPS.jl/wiki).

## Use

Refer to the instructions on this [wiki page](https://github.com/ProjectTorreyPines/SD4SOLPS.jl/wiki/Demo).

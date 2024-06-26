SHELL := /bin/zsh
help:
	@echo "Help Menu"
	@echo
	@echo "make env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories"
	@echo "make env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones"
	@echo "make clean: Deletes Project.toml and Manifest.toml for a fresh start"
	@echo

env_with_cloned_repo r:
	@echo "Pulling sample files using dvc"
	-dvc pull
	@echo "Creating Julia environment by creating local clones of dependent repositories"
	@echo "Cloning the repositories and generating Manifest.toml"
	-git clone "git@github.com:ProjectTorreyPines/IMASDD.jl.git" ../IMASDD; \
	git clone "git@github.com:ProjectTorreyPines/GGDUtils.jl.git" ../GGDUtils; \
	git clone "git@github.com:ProjectTorreyPines/SOLPS2IMAS.jl.git" ../SOLPS2IMAS; \
	julia --project=. -e 'using Pkg; Pkg.rm(["IMASDD", "GGDUtils", "SOLPS2IMAS", "EFIT"]); Pkg.develop(path="../IMASDD"); Pkg.develop(path="../GGDUtils"); Pkg.develop(path="../SOLPS2IMAS"); Pkg.add(url="git@github.com:JuliaFusion/EFIT.jl.git", rev="master"); Pkg.instantiate()'

env_with_git_url u:
	@echo "Pulling sample files using dvc"
	-dvc pull
	@echo "Creating Julia environment with the git urls without creating local clones"
	@echo "Generating Project.toml and Manifest.toml"
	julia --project=. -e 'using Pkg; Pkg.rm(["IMASDD", "GGDUtils", "SOLPS2IMAS", "EFIT"]); Pkg.add(url="git@github.com:ProjectTorreyPines/IMASDD.jl.git", rev="master"); Pkg.add(url="git@github.com:ProjectTorreyPines/GGDUtils.jl.git", rev="master"); Pkg.add(url="git@github.com:ProjectTorreyPines/SOLPS2IMAS.jl.git", rev="master"); Pkg.add(url="git@github.com:JuliaFusion/EFIT.jl.git", rev="master"); Pkg.instantiate()'

clean:
	@echo "Deleting Manifest.toml"
	- rm Manifest.toml

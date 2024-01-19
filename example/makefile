SHELL := /bin/zsh
help:
	@echo "Help Menu"
	@echo
	@echo "make env_with_cloned_repo: Creates a Julia environment with the cloned repositories"
	@echo "make env_with_git_url: Creates a Julia environment with the git urls without creating local clones"
	@echo "make clean: Deletes Project.toml and Manifest.toml for a fresh start"
	@echo "make Clean: Deletes Project.toml, Manifest.toml, and any cloned repositories for a fresh start"
	@echo

env_with_cloned_repo:
	chmod +x create_env_with_clone.sh
	./create_env_with_clone.sh

env_with_git_url:
	chmod +x create_env_with_git.sh
	./create_env_with_git.sh

clean:
	@echo "Deleting Project.toml and Manifest.toml"
	rm Project.toml Manifest.toml

Clean: clean
	chmod +x rm_cloned_repo.sh
	./rm_cloned_repo.sh
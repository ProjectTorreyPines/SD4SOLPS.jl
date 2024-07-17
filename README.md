# SD4SOLPS

![Format Check](https://github.com/ProjectTorreyPines/SD4SOLPS.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/SD4SOLPS.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/SD4SOLPS.jl/actions/workflows/test.yml/badge.svg)

Synthetic diagnostic workflow manager for use with SOLPS models

This repository is the top level layer for managing a workflow for calculating
the outputs of synthetic diagnostics attached to a SOLPS model.

Steps:
1) Load SOLPS into IMAS DD format if not already
2) Load equilibrium (that the SOLPS mesh was based on) into IMAS DD format
3) Make assumptions to extend profiles into the core and far SOL, if needed
4) Run synthetic diagnostic models and record output

For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/SD4SOLPS.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/SD4SOLPS.jl/dev).

## Installation

SD4SOLPS is registered with public repository [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). For installation:

```
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("SD4SOLPS")
```

## Examples

Refer to the instructions on this [wiki page](https://github.com/ProjectTorreyPines/SD4SOLPS.jl/wiki/Demo) to see how to run `examples/demo.ipynb`.

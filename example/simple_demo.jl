# A simple demo that can be run without jupyter notebooks.
# First, activate the project (preferably after loading Revise), then include this file:
# > using Revise
# ] activate .
# > include("example/simple_demo.jl")
using SOLPS2IMAS: SOLPS2IMAS
using GGDUtils: GGDUtils
using Plots
using SD4SOLPS

sample_path = "$(@__DIR__)/../sample/ITER_Lore_2296_00000/"
solps2imas_samples = splitdir(pathof(SOLPS2IMAS))[1] * "/../samples"

dd = SD4SOLPS.preparation(
    "Baseline2008-li0.70.x4.mod2.eqdsk",
    [sample_path, solps2imas_samples]...;
)

grid_ggd = dd.edge_profiles.grid_ggd[1] # First grid_ggd time slice. It is allowed to vary in time
space = grid_ggd.space[1] # First space in this grid_ggd

# Choose backend
gr()           # Fast and can save pdf
# plotlyjs()   # Use for interactive plot, can only save png
n_e = GGDUtils.get_prop_with_grid_subset_index(
    dd.edge_profiles.ggd[1].electrons.density,
    5,
)
plot(dd.edge_profiles.grid_ggd, n_e; colorbar_title="Electrons density / m^(-3)")
plot!(
    space,
    GGDUtils.get_grid_subset_with_index(grid_ggd, 16);
    linecolor=:black,
    linewidth=2,
    linestyle=:solid,
    label="Separatix",
    legend=true,
)

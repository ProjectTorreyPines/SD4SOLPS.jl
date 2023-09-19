"""
Utilities for extrapolating profiles
"""

# import CalculusWithJulia
import OMAS
import Interpolations
import NumericalIntegration
import GGDUtils

export extrapolate_core
export fill_in_extrapolated_core_profile
export mesh_psi_spacing

"""
    extrapolate_core(edge_rho, edge_quantity, rho_output)

Function for assuming a core profile when given edge profile data.

Concept:
1) value and derivative should be continuous when joining real data
2) derivative at magnetic axis is known to be 0 when making profile vs. rho,
    by the definition of rho
3) derivative probably does something fancier between the pedestal and axis than
    just linear interpolation, so add an extra point in there
4) there's a joint between the steep pedestal and the shallow core that needs an
    extra knot to manage it properly
5) after making up a very simple gradient profile out of a few line segments,
    integrate it to get the profile of the quantity in question
"""
function extrapolate_core(edge_rho, edge_quantity, rho_output)
    grad = OMAS.gradient(edge_rho, edge_quantity)
    gf = grad[1]
    rf = edge_rho[1]
    gmid = -abs(gf) / 4.0
    rmid = rf / 2.0
    rped_enforce = rf - 0.08
    gped_enforce = rf + (gmid - gf) / (rmid - rf) * rped_enforce
    gped_max = maximum(grad) / 10.0
    gped_enforce = minimum([abs(gped_enforce), abs(gped_max)]) * sign(gf)

    gg = [0, gmid, gped_enforce, gf]
    rr = [0, rmid, rped_enforce, rf]
    # https://www.matecdev.com/posts/julia-interpolation.html
    itp = Interpolations.LinearInterpolation(rr, gg, extrapolation_bc=Interpolations.Line())
    #itp = Interpolations.extrapolate(itp, Interpolations.Line())
    g = itp(rho_output)
    q_extend = NumericalIntegration.cumul_integrate(rho_output, g)
    q_offset = edge_quantity[1] - Interpolations.LinearInterpolation(rho_output, q_extend)(rf)
    q_extend .+= q_offset

    output_profile = Array{Float64}(undef, length(rho_output))
    output_profile[rho_output .< rf] = q_extend[rho_output .< rf]
    output_profile[rho_output .>= rf] = Interpolations.LinearInterpolation(edge_rho, edge_quantity).(rho_output[rho_output .>= rf])
    return output_profile
end

"""
    function fill_in_extrapolated_core_profile(dd::OMAS.dd, quantity_name::String)

This function accepts a DD that should be populated with equilibrium and edge_profiles
as well as a request for a quantity to extrapolate into the core. It then maps edge_profiles
data to rho, calls the function that performs the extrapolation (which is not a simple linear
extrapolation but has some trickery to attempt to make a somewhat convincing profile shape),
and writes the result to core_profiles. This involves a bunch of interpolations and stuff.

dd: an IMAS/OMAS data dictionary
quantity_name: the name of a quantity in edge_profiles.profiles_2d and core_profiles.profiles_1d,
    such as "electrons.density"
"""
function fill_in_extrapolated_core_profile(dd::OMAS.dd, quantity_name::String)
    println("hey hey hey")
    ggd_idx = 1
    space_number = 1
    space = dd.edge_profiles.grid_ggd[ggd_idx].space[space_number]
    cell_subset = SOLPS2IMAS.get_grid_subset_with_index(dd.edge_profiles.grid_ggd[ggd_idx], 5)
    midplane_subset = SOLPS2IMAS.get_grid_subset_with_index(dd.edge_profiles.grid_ggd[ggd_idx], 11)
    
    if length(midplane_subset.element) < 1
        throw(ArgumentError(string(
            "Midplane subset length was ", length(midplane_subset.element),
            ". Unacceptable data in data dictionary."
        )))
    end
    if length(cell_subset.element) < 1
        throw(ArgumentError(string(
            "Cell subset length was ", length(cell_subset.element),
            ". Unacceptable data in data dictionary."
        )))
    end

    nt = 1
    if length(dd.core_profiles.profiles_1d) < nt
        resize!(dd.core_profiles.profiles_1d, nt)
    end
    it = 1
    eq_it = 1
    # quantity_in_cells = SOLPS2IMAS.val_obj(dd.edge_profiles.ggd[it], quantity_name, ggd_idx)
    tags = split(quantity_name, ".")
    quantity_str = dd.edge_profiles.ggd[it]
    for tag in tags
        quantity_str = getproperty(quantity_str, Symbol(tag))
    end
    # Quantity is in cells now
    cell_subset_idx = 5
    nqv = length(quantity_str[cell_subset_idx].values)
    if nqv < 1
        println("Quantity ", quantity_name, " has ", nqv, " elements in subset ", cell_subset_idx, ". Unacceptable.")
        return
    end
    quantity_str[cell_subset_idx].grid_subset_index = cell_subset_idx
    midplane_cell_centers, quantity = GGDUtils.project_prop_on_subset!(quantity_str, cell_subset, midplane_subset, space)
    # Now quantity is at the outboard midplane

    # Get the rho values to go with the midplane quantity values
    # midplane_cell_centers = GGDUtils.get_subset_centers(space, midplane_subset)  # Not needed; done by project_prop_on_subset
    r = [midplane_cell_centers[i][1] for i in 1:length(midplane_cell_centers)]
    z = [midplane_cell_centers[i][2] for i in 1:length(midplane_cell_centers)]
    if length(r) != length(quantity)
        throw(DimensionMismatch(string(
            "Number of cell center coordinates (", length(r),
            ") does not match number of cell values (", length(quantity),
            ") in the midplane subset."
        )))
    end
    eq_idx = 1
    r_eq = dd.equilibrium.time_slice[eq_it].profiles_2d[eq_idx].grid.dim1
    z_eq = dd.equilibrium.time_slice[eq_it].profiles_2d[eq_idx].grid.dim2
    rho1_eq = dd.equilibrium.time_slice[eq_it].profiles_1d.rho_tor_norm
    psia = dd.equilibrium.time_slice[eq_it].global_quantities.psi_axis
    psib = dd.equilibrium.time_slice[eq_it].global_quantities.psi_boundary
    psi1_eq = (dd.equilibrium.time_slice[eq_it].profiles_1d.psi .- psia) ./ (psib - psia)
    psi2_eq = (dd.equilibrium.time_slice[eq_it].profiles_2d[eq_idx].psi .- psia) ./ (psib - psia)
    println(size(psi2_eq), ", ", size(r_eq), ", ", size(z_eq))
    rzpi = Interpolations.LinearInterpolation((z_eq, r_eq), psi2_eq)
    in_bounds = (r .< maximum(r_eq)) .& (r .> minimum(r_eq)) .& (z .> minimum(z_eq)) .& (z .< maximum(z_eq))
    println(in_bounds)
    psi_for_quantity = 10.0 .+ zeros(length(r))
    psi_for_quantity[in_bounds] = rzpi.(z[in_bounds], r[in_bounds])
    println(length(psi1_eq), ", ", length(rho1_eq))
    rho_for_quantity = copy(psi_for_quantity)
    in_bounds = psi_for_quantity .<= 1.0
    dpsi = diff(psi1_eq)
    drho = diff(rho1_eq)
    prepend!(dpsi, [0.0])
    prepend!(drho, [0.0])
    rho_for_quantity[in_bounds] = Interpolations.LinearInterpolation(psi1_eq, rho1_eq).(psi_for_quantity[in_bounds])
    
    # Make sure the output 1D rho grid exists; create it if needed
    if length(dd.core_profiles.profiles_1d[it].grid.rho_tor_norm) == 0
        resize!(dd.core_profiles.profiles_1d[it].grid.rho_tor_norm, 201)
        # If you don't like this default, then you should write grid.rho_tor_norm before calling this function.
        dd.core_profiles.profiles_1d[it].grid.rho_tor_norm = collect(LinRange(0, 1, 201))
    end
    rho_core = dd.core_profiles.profiles_1d[it].grid.rho_tor_norm

    # Finally, we're ready to call the extrapolation function and write the result
    quantity_core = extrapolate_core(rho_for_quantity, quantity, rho_core)
    parent = dd.core_profiles.profiles_1d[it]
    tags = split(quantity_name, ".")
    for tag in tags[1:end-1]
        parent = getproperty(parent, Symbol(tag))
    end
    setproperty!(parent, Symbol(tags[end]), quantity_core)
end

"""
    function extrapolate_edge_exp(quantity_edge::Vector{Float64}, dqdpsi::Vector{Float64}, psin_out::Vector{Float64})

Exponential decay version of edge profile extrapolation. Should work well for
many quantities, including Te and ne. If the exponential profile has no background
offset, then its amplitude and scale length can be completely defined by matching
the quantity value and first derivative at the edge of the mesh.

quantity_edge: values of some physics quantity in cells along the outer edge of the mesh
dqdpsi: Gradient of the quantity vs. psi, aligned perpendicular to the row of cells
    being used.
psin_out: Normalized psi values along a vector orthogonal to the row of cells along the edge.
    These psi_N values should be outside of the mesh (because the quantity is already known in the mesh)
    The output will be a matrix.
"""
function extrapolate_edge_exp(
    quantity_edge::Vector{Float64},
    dqdpsi::Vector{Float64},
    psin_out::Vector{Float64},
)
    x = psin_out - 1.0
    lambda = - quantity_edge / dqdpsi
    q0 = quantity_edge ./ exp(-x'./lambda)
    return y0 * exp(-x ./ lambda)
end

"""
    function mesh_psi_spacing(dd::OMAS.dd; eq_time_idx::Int64=1)

Inspects the mesh to see how far apart faces are in psi_N.
Requires that GGD and equilibrium are populated.

dd: a data dictionary instance with required data loaded into it
eq_time_idx: index of the equilibrium time slice to use. For a typical SOLPS run,
        the SOLPS mesh will be based on the equilibrium reconstruction at a single time,
        so the DD associated with the SOLPS run only needs one equilibrium time slice
        to be loaded. However, one could combine the complete equilibrium time series
        with the SOLPS run and then have to specify which slice of the equilibrium
        corresponds to the SOLPS mesh.
"""
function mesh_psi_spacing(dd::OMAS.dd; eq_time_idx::Int64=1)
    # Inspect input
    if length(dd.equilibrium.time_slice) < eq_time_idx
        throw(ArgumentError(string(
            "DD equilibrium does not have enough time slices: ",
            length(dd.equilibrium.time_slice), " slices were present but slice index ",
            eq_time_idx, " was requested.",
        )))
    end
    bad_ggd = (length(dd.edge_profiles.grid_ggd) < 1)
    if bad_ggd
        throw(ArgumentError(string(
            "Invalid GGD data."
        )))
    end

    # Get flux map
    eq_idx = 1  # Most cases won't have need for more than one of these, as far as I know
    eqt = dd.equilibrium.time_slice[eq_time_idx]
    p2 = eqt.profiles_2d[eq_idx]
    r_eq = p2.grid.dim1
    z_eq = p2.grid.dim2
    psi = p2.psi
    psia = eqt.global_quantities.psi_axis
    psib = eqt.global_quantities.psi_boundary
    psin_eq = (psi .- psia) ./ (psib - psia)
    rzpi = Interpolations.LinearInterpolation((z_eq, r_eq), psin_eq)
    println(minimum(r_eq), ", ", maximum(r_eq))
    println(minimum(z_eq), ", ", maximum(z_eq))

    # Get a row of cells. Since the mesh should be aligned to the flux surfaces,
    # it shouldn't matter which row is used, although the divertor rows might be
    # weird. So use the outboard midplane. That's always a solid choice.
    ggd_idx = 1
    space_number = 1
    space = dd.edge_profiles.grid_ggd[ggd_idx].space[space_number]
    midplane_subset = SOLPS2IMAS.get_grid_subset_with_index(dd.edge_profiles.grid_ggd[ggd_idx], 11)
    midplane_cell_centers = GGDUtils.get_subset_centers(space, midplane_subset)
    r_mesh = [midplane_cell_centers[i][1] for i in 1:length(midplane_cell_centers)]
    z_mesh = [midplane_cell_centers[i][2] for i in 1:length(midplane_cell_centers)]
    println(minimum(r_mesh), ", ", maximum(r_mesh))
    println(minimum(z_mesh), ", ", maximum(z_mesh))
    psin_mesh = rzpi.(z_mesh, r_mesh)
    # This should come out sorted, but with GGD, who knows.
    ii = sortperm(psin_mesh)
    psin = psin_mesh[ii]
    dpsin = psin[end] - psin[end-1]
    return dpsin
end

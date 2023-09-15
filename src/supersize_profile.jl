"""
Utilities for extrapolating profiles
"""

# import CalculusWithJulia
import OMAS
import Interpolations
import NumericalIntegration
import GGDUtils

function testytesty()
    println("yeah, the function ran")
    return nothing
end

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
    gped_max = maximum(grad) / 10.
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
    output_profile[rho_output .>= rf] = Interpolations.LinearInterpolation(edge_rho, edge_quantity)(rho_output[rho_output .>= rf])
    return output_profile
end

function fill_in_extrapolated_core_profile(dd::OMAS.dd, quantity_name::String)
    println("hey hey hey")
    ggd_idx = 1
    space_number = 1
    space = dd.edge_profiles.grid_ggd[ggd_idx].space[space_number]
    cell_subset = SOLPS2IMAS.get_grid_subset_with_index(dd.edge_profiles.grid_ggd[ggd_idx], 5)
    midplane_subset = SOLPS2IMAS.get_grid_subset_with_index(dd.edge_profiles.grid_ggd[ggd_idx], 11)
    
    if length(midplane_subset.element) < 1
        println("Midplane subset length was ", length(midplane_subset.element), ". Unacceptable.")
        return
    end
    if length(cell_subset.element) < 1
        println("Cell subset length was ", length(cell_subset.element), ". Unacceptable.")
        return
    end

    it = 1
    eq_it = 1
    # quantity_in_cells = SOLPS2IMAS.val_obj(dd.edge_profiles.ggd[it], quantity_name, ggd_idx)
    tags = split(quantity_name, ".")
    quantity = dd.edge_profiles.ggd[it]
    for tag in tags
        quantity = getproperty(quantity, Symbol(tag))
    end
    # Quantity is in cells now
    cell_subset_idx = 5
    nqv = length(quantity[cell_subset_idx].values)
    if nqv < 1
        println("Quantity ", quantity_name, " has ", nqv, " elements in subset ", cell_subset_idx, ". Unacceptable.")
        return
    end
    quantity[cell_subset_idx].grid_subset_index = cell_subset_idx
    GGDUtils.project_prop_on_subset!(quantity, cell_subset, midplane_subset, space=space)
    # Now quantity is at the outboard midplane

    midplane_cell_centers = GGDUtils.get_subset_centers(space, midplane_subset)
    r = [midplane_cell_centers[i][1] for i in 1:length(midplane_cell_centers)]
    z = [midplane_cell_centers[i][2] for i in 1:length(midplane_cell_centers)]
    println(midplane_cell_centers)
    eq_idx = 1
    r_eq = dd.equilibrium.time_slice[eq_it].profiles_2d[eq_idx].grid.dim1
    z_eq = dd.equilibrium.time_slice[eq_it].profiles_2d[eq_idx].grid.dim2
    rho1_eq = dd.equilibrium.time_slice[eq_it].profiles_1d.rho_tor_norm
    psi1_eq = dd.equilibrium.time_slice[eq_it].profiles_1d.psi
    psi2_eq = dd.equilibrium.time_slice[eq_it].profiles_2d[eq_idx].psi
    println(size(psi2_eq), ", ", size(r_eq), size(z_eq))
    psi_for_quantity = Interpolations.LinearInterpolation((z_eq, r_eq), psi2_eq).(z, r)
    println(length(psi1_eq), ", ", length(rho1_eq))
    rho_for_quantity = Interpolations.LinearInterpolation(psi1_eq, rho1_eq).(psi_for_quantity)
    
    rho_core = dd.core_profiles.profiles_1d[it].grid.rho_tor_norm
    quantity_core = Interpolations.LinearInterpolation(rho_for_quantity, quantity).(rho_core)
    parent = dd.core_profiles.profiles_1d[it]
    tags = split(quantity_name, ".")
    for tag in tags[1:end-1]
        parent = getproperty(parent, Symbol(tag))
    end
    setproperty!(parent, Symbol(tags[end]), quantity_core)
end
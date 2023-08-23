"""
Utilities for extrapolating profiles
"""

# import CalculusWithJulia
import OMAS
import Interpolations
import NumericalIntegration

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

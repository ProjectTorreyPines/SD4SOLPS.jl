using PhysicalConstants.CODATA2018: BoltzmannConstant, ElementaryCharge
using Unitful: Unitful

export gas_unit_converter

torrl_per_pam3 = Float64(
    Unitful.upreferred((1 * Unitful.Torr * Unitful.L) / (Unitful.Pa * Unitful.m^3)),
)
electrons_per_molecule = Dict(
    "H" => 1 * 2,  # Assumed to mean H2. How would you even puff a bunch of H1.
    "H2" => 1 * 2,
    "D" => 1 * 2,
    "D2" => 1 * 2,
    "T" => 1 * 2,
    "T2" => 1 * 2,
    "CH" => 6 + 1,
    "CD" => 6 + 1,
    "CH4" => 6 + 4,
    "CD4" => 6 + 4,
    "N" => 7 * 2,
    "N2" => 7 * 2,
    "Ne" => 10,
    "Ar" => 18,
    "Kr" => 36,
    "Xe" => 54,
)

"""
    gas_unit_converter(
        value_in::Float64,
        units_in::String,
        units_out::String;
        species::String="H",
        temperature::Float64=293.15,
    )

Converts gas flows between different units. Uses ideal gas law to convert between
Pressure * volume type flows / quantities and count / current types of units.
There is a version that accepts floats in and outputs floats, and another that
deals in Unitful quantities.
"""
function gas_unit_converter(
    value_in::Float64,
    units_in::String,
    units_out::String;
    species::String="H",
    temperature::Float64=293.15,
)
    if units_in == units_out
        return value_in
    end
    torrl_to_pam3 = torrl_per_pam3

    pam3_to_molecules = 1.0 / (temperature * BoltzmannConstant.val)
    torrl_to_molecules = torrl_to_pam3 * pam3_to_molecules

    factor_to_get_molecules_s = Dict(
        "torr L s^-1" => torrl_to_molecules,
        "molecules s^-1" => 1.0,
        "Pa m^3 s^-1" => pam3_to_molecules,
        "el s^-1" => 1.0 / electrons_per_molecule[species],
        "A" => ElementaryCharge / electrons_per_molecule[species],
    )
    factor_to_get_molecules = Dict(
        "torr L" => torrl_to_molecules,
        "molecules" => 1.0,
        "Pa m^3" => pam3_to_molecules,
        "el" => 1.0 / electrons_per_molecule[species],
        "C" => ElementaryCharge / electrons_per_molecule[species],
    )
    if haskey(factor_to_get_molecules_s, units_in) &
       haskey(factor_to_get_molecules_s, units_out)
        conversion_factor =
            factor_to_get_molecules_s[units_in] / factor_to_get_molecules_s[units_out]
    elseif haskey(factor_to_get_molecules, units_in) &
           haskey(factor_to_get_molecules, units_out)
        conversion_factor =
            factor_to_get_molecules[units_in] / factor_to_get_molecules[units_out]
    else
        throw(ArgumentError("Unrecognized units: " * units_in * " or " * units_out))
    end
    return value_in .* conversion_factor
end

"""
    gas_unit_converter(
        value_in::Unitful.Quantity,
        units_in::String,
        units_out::String;
        species::String="H",
        temperature=293.15 * Unitful.K,
    )

Converts gas flows between different units. Uses ideal gas law to convert between
Pressure * volume type flows / quantities and count / current types of units.
This is the Unitful version.

Output will be unitful, but the units are not simplified automatically. You can
perform operations such as

  - `(output |> Unitful.upreferred).val`
  - `Unitful.uconvert(Unitful.whatever, output).val`

to handle simplification or conversion of units.

Although this function pretends torr L s``^{-1}`` and Pa m``^3`` s``^{-1}`` are
different, use of Unitful should cause them to behave the same way as long as you
simplify or convert units at the end. This means that you can use other pressure*volume
type gas units and call them torr L s``^{-1}`` and the script will deal with them up to
having messy units in the output.
"""
function gas_unit_converter(
    value_in::Unitful.Quantity,
    units_in::String,
    units_out::String;
    species::String="H",
    temperature=293.15 * Unitful.K,
)
    if units_in == units_out
        return value_in
    end

    # The conversion between Torr*L and Pa*m^-3 is the only mundane unit conversion in
    # here. Unitful or any other unit conversion tool could handle this by itself
    # (e.g. Unitful.uconvert), but you wouldn't need this fancy function for that.
    # The rest of the conversions, the hard ones (that require assuming a temperature
    # or knowing molecule structure and atomic numbers) are all implemented with
    # multiplicative factors, so Unitful's normal conversion handling is bypassed
    # and replaced with a multiplicative factor.
    torrl_to_pam3 = torrl_per_pam3 * Unitful.Pa * Unitful.m^3 / Unitful.Torr / Unitful.L
    # torrl_to_pam3 should simplify to 1, but Unitful doesn't do this simpliciation
    # without prompting. So the Torr L or the Pa m^-3 will cancel when this factor
    # is multiplied by some gas pressure*volume units.

    pam3_to_molecules =
        Unitful.J / (temperature * BoltzmannConstant) / (Unitful.Pa * Unitful.m^3)
    torrl_to_molecules = torrl_to_pam3 * pam3_to_molecules

    factor_to_get_molecules_s = Dict(
        "torr L s^-1" => torrl_to_molecules,
        "molecules s^-1" => 1.0,
        "Pa m^3 s^-1" => pam3_to_molecules,
        "el s^-1" => 1.0 / electrons_per_molecule[species],
        "A" => ElementaryCharge / electrons_per_molecule[species],
    )
    factor_to_get_molecules = Dict(
        "torr L" => torrl_to_molecules,
        "molecules" => 1.0,
        "Pa m^3" => pam3_to_molecules,
        "el" => 1.0 / electrons_per_molecule[species],
        "C" => ElementaryCharge / electrons_per_molecule[species],
    )
    if haskey(factor_to_get_molecules_s, units_in) &
       haskey(factor_to_get_molecules_s, units_out)
        conversion_factor =
            factor_to_get_molecules_s[units_in] / factor_to_get_molecules_s[units_out]
    elseif haskey(factor_to_get_molecules, units_in) &
           haskey(factor_to_get_molecules, units_out)
        conversion_factor =
            factor_to_get_molecules[units_in] / factor_to_get_molecules[units_out]
    else
        throw(ArgumentError("Unrecognized units: " * units_in * " or " * units_out))
    end
    return value_in .* conversion_factor
end

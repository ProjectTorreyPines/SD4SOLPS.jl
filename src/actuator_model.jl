# Actuator models to translate commands (probably in V) into gas flows

using PhysicalConstants.CODATA2018
using Unitful: Unitful
using Interpolations: Interpolations
using YAML: YAML

"""
gas_unit_converter()

Converts gas flows between different units. Uses ideal gas law to convert between
Pressure * volume type flows / quantities and count / current types of units.
Output will be unitful; take output.val if you do no want the units attached.
"""
function gas_unit_converter(
    value_in, units_in, units_out; species::String="H", temperature=293.15 * Unitful.K,
)
    if units_in == units_out
        return value_in
    end

    # Multiply by gas flow to convert Torr L/s to Pa m^3/s
    torrl_to_pam3 = 0.133322368 * Unitful.Pa * Unitful.m^3 / Unitful.Torr / Unitful.L
    pam3_to_molecules =
        Unitful.J / (temperature * BoltzmannConstant) / (Unitful.Pa * Unitful.m^3)
    torrl_to_molecules = torrl_to_pam3 * pam3_to_molecules
    atoms_per_molecule = Dict(
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
    atoms_per_molecule[species]

    factor_to_get_molecules_s = Dict(
        "torr L s^-1" => torrl_to_molecules,
        "molecules s^-1" => 1.0,
        "Pa m^3 s^-1" => pam3_to_molecules,
        "el s^-1" => 1.0 / atoms_per_molecule[species],
        "A" => ElementaryCharge / atoms_per_molecule[species],
    )
    factor_to_get_molecules = Dict(
        "torr L" => torrl_to_molecules,
        "molecules" => 1.0,
        "Pa m^3" => pam3_to_molecules,
        "el" => 1.0 / atoms_per_molecule[species],
        "C" => ElementaryCharge / atoms_per_molecule[species],
    )
    if haskey(factor_to_get_molecules_s, units_in)
        conversion_factor =
            factor_to_get_molecules_s[units_in] / factor_to_get_molecules_s[units_out]
    elseif haskey(factor_to_get_molecules, units_in)
        conversion_factor =
            factor_to_get_molecules[units_in] / factor_to_get_molecules[units_out]
    else
        throw(ArgumentError("Unrecognized units: " * units_in))
    end
    return value_in .* conversion_factor
end

select_default_config(model::String) = model * "_gas_valve.yml"

"""
    model_gas_valve()

The main function for handling a gas valve model. Has logic for selecting models
and configurations.
"""
function model_gas_valve(
    t, command, model::String; configuration_file::String="auto", species::String="D2",
)
    # Select configuration file
    if configuration_file == "auto"
        configuration_file = select_default_config(model)
    end
    default_config_path = "$(@__DIR__)/../config/"
    if !occursin("/", configuration_file)
        configuration_file = default_config_path * configuration_file
    end
    if !isfile(configuration_file)
        throw(ArgumentError("Configuration file not found: " * configuration_file))
    end
    config = YAML.load_file(configuration_file)
    if (species in ["H", "H2", "D", "T", "T2"]) & !(species in keys(config))
        config = config["D2"]  # They should all be the same
    elseif (species in ["CH"]) & !(species in keys(config))
        config = config["CD"]
    elseif (species in ["CH4", "CD4"]) & !(species in keys(config))
        config = config["CD4"]
    else
        config = config[species]
    end

    if model == "simple"
        return simple_gas_model(t, command, config)
    else
        throw(ArgumentError("Unrecognized model: " * model))
    end
end

function instant_gas_model(command, config)
    return config["p1"] .* (sqrt.(((command * config["p2"]) .^ 2.0 .+ 1) .- 1))
end

function lowpass_filter_(raw, previous_smooth, dt, tau)
    return previous_smooth + dt / (dt + tau) * (raw - previous_smooth)
end

function lowpass_filter(t, x, tau)
    xs = zeros(length(t))
    for i ∈ 2:length(t)
        xs[i] = lowpass_filter_(x[i], xs[i-1], t[i] - t[i-1], tau)
    end
    return xs
end

function simple_gas_model(t, command, config)
    flow0 = instant_gas_model(command, config)
    t_ext = copy(t)
    prepend!(t_ext, minimum(t) - config["delay"])
    flow0_ext = copy(flow0)
    prepend!(flow0_ext, flow0[1])
    interp = Interpolations.LinearInterpolation(t_ext, flow0_ext)
    delayed_flow = interp.(t .- config["delay"])
    return flow = lowpass_filter(t, delayed_flow, config["tau"])
end

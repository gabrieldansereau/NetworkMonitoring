using DrWatson
@quickactivate :NetworkMonitoring

# Random.seed!(42)

# Define all required variables
# set_params = true
const defaults = Dict(
    :ns => 75,
    :nsites => 100,
    :C_exp => 0.2,
    :ra_sigma => 1.2,
    :ra_scaling => 50.0,
    :energy_NFL => 10_000,
    :H_nlm => 0.5,
    :nbon => 50,
    :refmethod => "metawebify",
)

## Generate network & abundances

function generate_networks(;
    ns=defaults[:ns],
    nsites=defaults[:nsites],
    C_exp=defaults[:C_exp],
    ra_sigma=defaults[:ra_sigma],
    ra_scaling=defaults[:ra_scaling],
    energy_NFL=defaults[:energy_NFL],
)
    # Generate metaweb using Niche Model
    # Random.seed!(42);
    metaweb = generate(SIS.NicheModel(ns, C_exp))

    # Generate autocorrelated ranges
    ranges = generate(AutocorrelatedRange(dims=(nsites, nsites)),ns)

    # Generate realized abundances
    ra = generate(NormalizedLogNormal(σ=ra_sigma), metaweb)

    # Generate possible local networks
    pos = possible(metaweb, ranges)

    # Generate detectable metaweb given abundances
    detectable = detectability(RelativeAbundanceScaled(ra_scaling), metaweb, ra)

    # Generate possible spatial metaweb (given cooccurrence)
    realized = possible(metaweb, ranges)
    # Make network realizable given Neutrally forbidden links
    realizable!(NeutrallyForbiddenLinks(energy_NFL), realized, ra)
    # Realize it
    realize!(realized)

    # Detect it
    detected = deepcopy(realized)
    detect!(detected, detectable)

    res = @dict metaweb ranges ra pos realized detected

    return res
end

## Add BON

function generate_bon(;
    nsites=defaults[:nsites],
    H_nlm=defaults[:H_nlm],
    nbon=defaults[:nbon],
)
    # Generate uncertainty layer
    uncertainty = SDT.SDMLayer(
        MidpointDisplacement(H_nlm), (nsites, nsites);
        x=(0.0, nsites), y=(0.0, nsites)
    )

    # Use BalancedAcceptance sampling
    bon = BON.sample(BON.BalancedAcceptance(nbon), uncertainty)

    return bon
end

## Extract infos from BON & evaluate monitoring

# Species
function evaluate_monitoring(
    occ::T, bon::BON.BiodiversityObservationNetwork
) where T <: Occurrence
    # List monitored species
    monitored = monitor(x -> findall(isone, x), occ, bon; makeunique=true)
    monitored = NetworkMonitoring._getspecies(monitored, occ)

    # Evaluate proportion of monitored species
    prop_monitored = length(monitored) / length(occ.species)
end

# Interactions
function evaluate_monitoring(
    m::T, bon::BON.BiodiversityObservationNetwork; ref=nothing
) where T <: Metaweb
    # List monitored interactions
    monitored = monitor(x -> interactions(render(Binary, x)), m, bon; makeunique=true)

    # Rebuild global metaweb for reference
    m_global = metawebify(m)

    # List proportions to compare
    if isnothing(ref)
        prop_monitored = length(monitored) / sum(m_global)
    else
        prop_monitored = length(monitored) / ref
    end
end

## Run all

function main(d::Dict)
    # Extract parameters
    @unpack ns, nsites, C_exp, ra_sigma, ra_scaling, energy_NFL, H_nlm, nbon, refmethod = d

    # Generate networks using simulations
    nets_dict = generate_networks(;
        ns=ns,
        nsites=nsites,
        C_exp=C_exp,
        ra_sigma=ra_sigma,
        ra_scaling=ra_scaling,
        energy_NFL=energy_NFL,
    )
    @unpack detected, realized, pos, metaweb, ranges = nets_dict

    # Generate BON
    bon = generate_bon(;
        nsites=nsites,
        H_nlm=H_nlm,
        nbon=nbon,
    )

    # Evaluate species monitoring
    prop_monitored_sp = evaluate_monitoring(ranges, bon)

    # Evaluate interactions monitoring
    if refmethod == "metawebify"
        ref = nothing
    elseif refmethod == "global"
        ref = length(metaweb.metaweb)
    end
    prop_detected_int = evaluate_monitoring(detected, bon; ref=ref)
    prop_realized_int = evaluate_monitoring(realized, bon; ref=ref)
    prop_possible_int = evaluate_monitoring(pos, bon; ref=ref)

    return @dict prop_detected_int prop_realized_int prop_possible_int prop_monitored_sp
end
main() = main(defaults)

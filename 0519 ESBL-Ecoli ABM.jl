#ESBL-producing E. coli Dynamic ABM in Tanzanian NICUs (version 1.0)
# Luyuan Zhao - 2025-05-19 

# 1. Loading necesaary packages
using Random               # random number generation
using StatsBase            # sample without replacement
using Plots;
gr();                      # plotting backend

# 2. State & agent enumerations
@enum BabyState S C_S C_R I_S I_R REC
# To track the infection or AMR status of each neonate
# S=Susceptible, C_S=Colonized Susceptible, C_R=Colonized Resistant, 
# I_S=Infected Susceptible, I_R=Infected Resistant, REC=Recovered

@enum CarrierState CLEAN CONTAM
# Describe whether a health-care worker (HCW) or device is contaminated

@enum AgentType BABY HCW DEVICE
# Let one common struct represent babies, HCWs and devices while still knowing “what kind of thing am I?”

# 3. Agent & Environment settings
# Agent - Stores identification, type, health status, contamination status and timestamps for epidemiological events.
mutable struct Agent
    id::Int
    atype::AgentType
    bstate::BabyState      # only meaningful for babies
    cstate::CarrierState   # only meaningful for HCW / device
    variant::Symbol        # :S or :R
    tick::Int              # days since entering current state
    abx_left::Int          # antibiotic days remaining
end

mutable struct Environment
    babies::Vector{Agent}
    hcws::Vector{Agent}
    devices::Vector{Agent}

    # All βeta parameters are the probabilities of transmission between agents.
    β_hb::Float64          # HCW -> Baby
    β_bh::Float64          # Baby -> HCW
    β_db::Float64          # Device -> Baby
    β_bd::Float64          # Baby -> Device
    β_dh::Float64          # Device -> HCW
    β_hd::Float64          # HCW -> Device  

    # The number of contacts among babies, HCWs and devices.
    contacts_hb::Int       # HCW <-> Baby contacts per day
    contacts_db::Int       # Device <-> Baby contacts
    contacts_dh::Int       # HCW <-> Device contacts  

    # The probability of all bacteria changes in the environment.
    p_progress::Float64    # colonised -> infection
    p_recovery_S::Float64  # I_S recovery
    p_recovery_R::Float64  # I_R recovery
    abx_course::Int        # length of empiric therapy (days)
    p_selection::Float64   # probability C_S -> C_R under ABX

    # The probability of decontamination of HCWs and devices.
    p_decontam_hcw::Float64
    p_decontam_dev::Float64

    # Time steps, all stats needed 
    # Counts the total number of antibiotic treatment days given to babies over the whole simulation
    # and the number of selection events (i.e. when a baby is given antibiotics).
    tick::Int
    stats::Dict{BabyState,Vector{Int}}
    total_abx_days::Int
    selection_events::Int
end

# 4.  Initialization
# Initial Status (default values)
# There are 10 babies, 3 HCWs and 5 devices in the environment.
function init_env(; N_babies=10, N_hcws=3, N_devices=5,
    β_hb=0.04, β_bh=0.12, β_db=0.08, β_bd=0.2,
    β_dh=0.05, β_hd=0.05,
    contacts_hb=67, contacts_db=20, contacts_dh=30,
    p_progress=0.03, p_recovery_S=0.20, p_recovery_R=0.12,
    abx_course=5, p_selection=0.25,
    p_decontam_hcw=0.40, p_decontam_dev=0.30,
    rng=Random.GLOBAL_RNG)

    # We introduce 1 initial contamination or colonisation event in each agent type
    #(baby, healthcare worker, and medical device) to avoid a sterile system.
    # Create babies
    babies = Agent[]
    for i in 1:N_babies
        bstate = (i == 1) ? C_S : S          # first baby starts colonized
        baby = Agent(i, BABY, bstate, CLEAN, :S, 0, 0)
        push!(babies, baby)
    end

    # Create HCWs
    hcws = Agent[]
    for i in 1:N_hcws
        hcw = Agent(i, HCW, S, CLEAN, :S, 0, 0)
        push!(hcws, hcw)
    end

    # Create devices
    devices = Agent[]
    for i in 1:N_devices
        device = Agent(i, DEVICE, S, CLEAN, :S, 0, 0)
        push!(devices, device)
    end

    # stats dict initialisation
    stats = Dict(s => Int[] for s in instances(BabyState))
    

    return Environment(babies, hcws, devices,
        β_hb, β_bh, β_db, β_bd, β_dh, β_hd,
        contacts_hb, contacts_db, contacts_dh,
        p_progress, p_recovery_S, p_recovery_R,
        abx_course, p_selection,
        p_decontam_hcw, p_decontam_dev,
        0, stats, 0, 0)
end

# 5. Helper: sample without replacement
sample_norpl(rng, v::Vector, n::Int) = StatsBase.sample(rng, v, min(n, length(v)), replace=false)

# 6. Record stats
# To record the number of agents in each state and add it to the stats dictionary.
function record_stats!(env::Environment)
    counts = countmap(b.bstate for b in env.babies)
    for s in instances(BabyState)
        push!(env.stats[s], get(counts, s, 0))
    end
end

# 7. Daily dynamics
# The function daily_step!() simulates the daily dynamics of the environment.
function daily_step!(env::Environment, rng::AbstractRNG)

    #HCW <=> Baby contacts (bidirectional) ----------
    for h in env.hcws
        for b in sample_norpl(rng, env.babies, env.contacts_hb)

            # HCW -> Baby: A contaminated HCW touches a susceptible baby
            # with transmission probability β_hb, baby becomes colonized with susceptible strain
            if h.cstate == CONTAM && b.bstate == S && rand(rng) < env.β_hb
                b.bstate = C_S
                b.variant = :S
                b.tick = 0
            end

            # Baby -> HCW: A colonized/infected baby contaminates a clean HCW
            # with probability β_bh during contact
            if b.bstate in (C_S, C_R, I_S, I_R) && h.cstate == CLEAN && rand(rng) < env.β_bh
                h.cstate = CONTAM
            end
        end
    end

    # Device ↔ Baby contacts
    for d in env.devices
        for b in sample_norpl(rng, env.babies, env.contacts_db)

            # Device -> Baby: A contaminated device touches a susceptible baby
            # with transmission probability β_db
            if d.cstate == CONTAM && b.bstate == S && rand(rng) < env.β_db
                b.bstate = C_S
                b.variant = :S
                b.tick = 0
            end

            # Baby -> Device: A colonized/infected baby contaminates a clean device
            # with probability β_bd during contact
            if b.bstate in (C_S, C_R, I_S, I_R) && d.cstate == CLEAN && rand(rng) < env.β_bd
                d.cstate = CONTAM
            end
        end
    end

    # HCW <=> Device contacts (NEW bidirectional) ---
    for h in env.hcws
        for d in sample_norpl(rng, env.devices, env.contacts_dh)
            # Device → HCW: A contaminated device contaminates a clean HCW
            if d.cstate == CONTAM && h.cstate == CLEAN && rand(rng) < env.β_dh
                h.cstate = CONTAM
            end

            # HCW -> Device: A contaminated HCW contaminates a clean device 
            if h.cstate == CONTAM && d.cstate == CLEAN && rand(rng) < env.β_hd
                d.cstate = CONTAM
            end
        end
    end

    # HCW / Device decontamination 
    # Each contaminated HCW has a chance to become clean (hand hygiene)
    for h in env.hcws
        h.cstate = (h.cstate == CONTAM && rand(rng) < env.p_decontam_hcw) ? CLEAN : h.cstate
    end
    # Each contaminated device has a chance to be disinfected
    for d in env.devices
        d.cstate = (d.cstate == CONTAM && rand(rng) < env.p_decontam_dev) ? CLEAN : d.cstate
    end

    # Within-baby events 
    for b in env.babies

        # Antibiotic usage: if currently receiving antibiotics
        if b.abx_left > 0
            b.abx_left -= 1
            env.total_abx_days += 1

            # Selection pressure: C_S -> C_R during antibiotic treatment
            if b.bstate == C_S && rand(rng) < env.p_selection
                b.bstate = C_R
                b.variant = :R
                b.tick = 0
                env.selection_events += 1
            end
        end

        # Disease progression: colonization becomes infection
        if b.bstate in (C_S, C_R) && rand(rng) < env.p_progress
            b.bstate = (b.bstate == C_S) ? I_S : I_R
            b.abx_left = env.abx_course
            b.tick = 0
        end

        # Recovery from infection: I_S or I_R → REC
        if (b.bstate == I_S && rand(rng) < env.p_recovery_S) ||
           (b.bstate == I_R && rand(rng) < env.p_recovery_R)
            b.bstate = REC
            b.abx_left = 0
            b.tick = 0
        end

        b.tick += 1
    end

    env.tick += 1
end

# 8. Simulation Loop
# Simulation controller
# Runs the full simulation for a given number of days (90 days here by default)
# and record the number of babies in each state (S, C_S, C_R, I_S, I_R, REC)
function run!(env::Environment; days=90, rng=Random.GLOBAL_RNG)
    record_stats!(env)
    for _ in 1:days
        daily_step!(env, rng)
        record_stats!(env)
    end
end

# 9. Results and Visualisations
function main()
    println("Starting E. coli ABM (with initial contamination)...")

    env = init_env()
    run!(env, days=90)

    # Compute total infections per day
    infections = [env.stats[I_S][t] + env.stats[I_R][t] for t in 1:length(env.stats[I_S])]
    peak_inf = maximum(infections)

    # Plotting
    days = 0:(length(env.stats[S])-1)
    plt = plot(title="ESBL-producing E. coli Dynamics in NICU (by default)",
        xlabel="Days", ylabel="Neonate Counts", legend=:outertopright)

    for s in instances(BabyState)
        plot!(plt, days, env.stats[s], label=string(s))
    end
    display(plt)

    # Summary output
    println("\nPeak infections (I_S + I_R): ", peak_inf)
    println("Total antibiotic days: ", env.total_abx_days)
    println("C_S → C_R selection events:  ", env.selection_events)

    # Average ABX days per baby
    # (total ABX days / number of babies)
    N_babies = length(env.babies)
    avg_abx_days_per_baby = round(env.total_abx_days / N_babies, digits=2)
    println("Avg total ABX treatment days per baby: ", avg_abx_days_per_baby)
end

main()
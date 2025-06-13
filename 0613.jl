# ESBL-producing Enterobacterales Dynamic ABM in Tanzanian NICUs (version 1.7)
# Luyuan Zhao 
# 2025-06-13 
# In this version, we added the expected outputs in the baseline model.


#---------------------------------------------------------------------------------------------------------------------
# 1. Loading necessary packages
using Random               # random number generation
using StatsBase            # sample without replacement
using Plots;               # for visualization
using Statistics           # for mean and ± SD    
using Measures             # for units
gr();                      # plotting backend


#---------------------------------------------------------------------------------------------------------------------
  # 2. State & agent enumerations
@enum BabyState S C_S C_R I_S I_R REC
  # To track the infection or AMR status of each neonate: 
  #  S = Susceptible, C_S = Colonized Susceptible, C_R = Colonized Resistant, 
  #  I_S = Infected Susceptible, I_R = Infected Resistant, REC = Recovered

@enum CarrierState CLEAN CONTAM
  # Describe whether a health-care worker (HCW) or device is contaminated

@enum AgentType BABY HCW DEVICE
  # Let one common struct represent babies, HCWs and devices while still knowing “what kind of thing am I?”


#---------------------------------------------------------------------------------------------------------------------
# 3. Agent & Environment settings
# 3.1 Agent - Stores identification, type, health status, contamination status and timestamps for epidemiological events.
mutable struct Agent
    id::Int
    atype::AgentType
    bstate::BabyState      # only meaningful for babies
    cstate::CarrierState   # only meaningful for HCW / device
    variant::Symbol        # :S or :R
    tick::Int              # days since entering current state
    abx_left::Int          # antibiotic days remaining
end

# 3.2 Environment (Mutable)- Stores all agents and parameters of the simulation.
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
    p_progress::Float64    # colonisation -> infection
    p_recovery_S::Float64  # recovery probability (I_S -> R)
    p_recovery_R::Float64  # recovery probability (I_R -> R)
    abx_course::Int        # length of empiric therapy (days)
    p_selection::Float64   # probability of selection pressure (C_S -> C_R) during antibiotic treatment

    # The probability of decontamination of HCWs and devices.
    p_decontam_hcw::Float64 # HCW can become clean
    p_decontam_dev::Float64 # Devices can become clean

    tick::Int                             # Time step
    stats::Dict{BabyState,Vector{Int}}    # Stats dictionary
    total_abx_days::Int                   # Total antibiotic days (i.e. days of treatment)
    selection_events::Int                 # Number of selection events
    total_recoveries::Int                 # Total number of recoveries

    p_empiric_abx::Float64                # probability of offering empiric therapy to babies

    p_exit::Float64               # Probability of leaving the NICU (discharge or death) per day
    p_init_C_S::Float64           # Probability of new babies being colonised (sensitive bacteria)
    p_init_C_R::Float64           # ~ being colonised (resistant bacteria)
    p_init_I_S::Float64           # Probability of new babies being infected (sensitive bacteria)
    p_init_I_R::Float64           # ~ being infected (resistant bacteria)
    baby_pool::Vector{Int}        # Pool of IDs for new babies
end


#---------------------------------------------------------------------------------------------------------------------
# 4.  Initialization (giving initial status or default values)
# There are always 10 babies, 3 HCWs and 5 devices in the environment.
function init_env(; N_babies=10, N_total=165, N_hcws=3, N_devices=5,
    β_hb=0.04, β_bh=0.12, β_db=0.08, β_bd=0.2,
    β_dh=0.05, β_hd=0.05,
    contacts_hb=67, contacts_db=20, contacts_dh=30,
    p_progress=0.03, p_recovery_S=0.20, p_recovery_R=0.12,
    abx_course=5, p_selection=0.25,
    p_decontam_hcw=0.40, p_decontam_dev=0.30, p_empiric_abx=0.15,
    p_exit=0.12,
    p_init_C_S=0.05, p_init_C_R=0.0,
    p_init_I_S=0.02, p_init_I_R=0.006,
    _=Random.GLOBAL_RNG)

 # Build up pools for babies
 all_ids = collect(1:N_total)
 init_ids = all_ids[1:N_babies]           # the first 10 babies are initialised with specific states
 baby_pool = all_ids[(N_babies+1):end]    # the rest of the IDs are available for new babies

    #  We introduce 3 initial cases: 1 I_R, 1 I_S and 1 C_S.
    # Create babies 
    babies = Agent[]
    for i in init_ids
        bstate, abx_left =
            if i == 1
                I_R, round(Int, abx_course * 1.5)    # I_R would receive 1.5x treatment duration
            elseif i == 2
                I_S, abx_course                      # I_S would have normal treatment
            elseif i == 3
                C_S, 0                               # Theoretically, C_S would not receive treatment
            else
                S, 0
            end
        baby = Agent(i, BABY, bstate, CLEAN, :S, 0, abx_left)  # id, type, state, contamination status, variant, tick, abx_left
        push!(babies, baby)
    end

    # Create HCWs 
    hcws = Agent[]
    for i in 1:N_hcws
        hcw = Agent(i, HCW, S, CLEAN, :S, 0, 0)         # id, type, state, contamination status, variant, tick, abx_left 
        push!(hcws, hcw)
    end

    # Create devices
    devices = Agent[]
    for i in 1:N_devices
        device = Agent(i, DEVICE, S, CLEAN, :S, 0, 0)   # id, type, state, contamination status, variant, tick, abx_left
        push!(devices, device)
    end

    # stats dict initialisation
    stats = Dict(s => Int[] for s in instances(BabyState))  # Create a dictionary to store the stats of each state

    return Environment(babies, hcws, devices,
        β_hb, β_bh, β_db, β_bd, β_dh, β_hd,
        contacts_hb, contacts_db, contacts_dh,
        p_progress, p_recovery_S, p_recovery_R,
        abx_course, p_selection,
        p_decontam_hcw, p_decontam_dev,
        0, stats, 0, 0, 0,
        p_empiric_abx, p_exit::Float64,              
        p_init_C_S::Float64,          
        p_init_C_R::Float64,          
        p_init_I_S::Float64,          
        p_init_I_R::Float64, 
        baby_pool)
end


#---------------------------------------------------------------------------------------------------------------------
# 5. Helper: sample without replacement
sample_norpl(rng, v::Vector, n::Int) = StatsBase.sample(rng, v, min(n, length(v)), replace=false)


#---------------------------------------------------------------------------------------------------------------------
# 6. Record stats (to record the number of agents in each state and add it to the stats dictionary)
function record_stats!(env::Environment)
    # Count the number of babies in each state and add it to the stats dictionary.
    counts = countmap(b.bstate for b in env.babies)

    # Record all possible status (6) traversed
    for s in instances(BabyState)

        # Add the number of people in this status today to the corresponding time series.
        push!(env.stats[s], get(counts, s, 0))
    end

end


#---------------------------------------------------------------------------------------------------------------------
# 7. Daily dynamics (simulate how this AMR transmission occurs in the NICU)
function daily_step!(env::Environment, rng::AbstractRNG)

    # HCW <=> Baby contacts
    for h in env.hcws
        for b in sample_norpl(rng, env.babies, env.contacts_hb)

            # HCW -> Baby: A contaminated HCW touches a susceptible baby
            # with transmission probability β_hb, baby becomes C_S or C_R
            if h.cstate == CONTAM && b.bstate == S && rand(rng) < env.β_hb
                if h.variant == :S
                    b.bstate = C_S
                elseif h.variant == :R
                    b.bstate = C_R
                end
                b.variant = h.variant
                b.tick = 0
            end

            # Baby -> HCW: A contaminated baby contaminates a clean HCW
            # with probability β_bh during contact
            if b.bstate in (C_S, C_R, I_S, I_R) && h.cstate == CLEAN && rand(rng) < env.β_bh
                h.variant = b.variant
            end
        end
    end

    # Device <=> Baby contacts
    for d in env.devices
        for b in sample_norpl(rng, env.babies, env.contacts_db)

            # Device -> Baby: a contaminated Device touches a susceptible baby
            # with transmission probability β_db, baby becomes C_S or C_R
            if d.cstate == CONTAM && b.bstate == S && rand(rng) < env.β_db
                if d.variant == :S
                    b.bstate = C_S
                elseif d.variant == :R
                    b.bstate = C_R
                end
                b.variant = d.variant
                b.tick = 0
            end

            # Baby -> Device: a contaminated baby contaminates a clean device
            # with probability β_bd during contact
            if b.bstate in (C_S, C_R, I_S, I_R) && d.cstate == CLEAN && rand(rng) < env.β_bd
                d.cstate = CONTAM
                d.variant = b.variant
            end
        end
    end

    # HCW <=> Device contacts (bidirectional)
    for h in env.hcws
        for d in sample_norpl(rng, env.devices, env.contacts_dh)
            # HCW -> Device
            if h.cstate == CONTAM && d.cstate == CLEAN && rand(rng) < env.β_hd
                d.cstate = CONTAM
                d.variant = h.variant
            end

            # Device -> HCW
            if d.cstate == CONTAM && h.cstate == CLEAN && rand(rng) < env.β_dh
                h.cstate = CONTAM
                h.variant = d.variant
            end
        end
    end

    # HCW or Device decontamination 
    # Each contaminated HCW has a chance to become clean (hand hygiene)
    for h in env.hcws
        if h.cstate == CONTAM && rand(rng) < env.p_decontam_hcw
            h.cstate = CLEAN
            h.variant = :S
        end
    end

    # Each contaminated device has a chance to be disinfected
    for d in env.devices
        if d.cstate == CONTAM && rand(rng) < env.p_decontam_dev
            d.cstate = CLEAN
            d.variant = :S
        end
    end

    # Within-baby events 
    for b in env.babies

        # 1) empirical ABX for colonised susceptible babies (prevention purpose)
        if b.bstate == C_S && b.abx_left == 0 && rand(rng) < env.p_empiric_abx
            b.abx_left = env.abx_course
        end

        # 2) apply antibiotics & selection pressure
        if b.abx_left > 0
            b.abx_left -= 1
            env.total_abx_days += 1
            if b.bstate == C_S && rand(rng) < env.p_selection
                b.bstate = C_R
                b.variant = :R
                b.tick = 0
                env.selection_events += 1
            end
        end

        # 3) progression to infection
        if b.bstate == C_S && rand(rng) < env.p_progress
            b.bstate = I_S
            b.abx_left = env.abx_course
            b.tick = 0
        elseif b.bstate == C_R && rand(rng) < env.p_progress
            b.bstate = I_R
            b.abx_left = env.abx_course
            b.tick = 0
        end

        # 4) recovery
        if (b.bstate == I_S && rand(rng) < env.p_recovery_S) ||
           (b.bstate == I_R && rand(rng) < env.p_recovery_R)
            b.bstate = REC
            b.abx_left = 0
            b.tick = 0

            env.total_recoveries += 1
            # REC  ->  S
            if b.bstate == REC
                b.bstate = S
            end
        end

        b.tick += 1
    end

    # Rules for babies leaving the NICU
    # find out which babies are discharged
    discharged = [b for b in env.babies if rand(rng) < env.p_exit]

    # rid them of the environment
    for b in discharged
        deleteat!(env.babies, findfirst(x -> x.id == b.id, env.babies))
    end

    # new addmissions
    # If there are less than 10 babies, we add new ones from the baby pool. 
    while length(env.babies) < 10 && !isempty(env.baby_pool)
        new_id = popfirst!(env.baby_pool)

        r = rand(rng)
        if r < env.p_init_I_R
            bstate, variant, abx = I_R, :R, round(Int, env.abx_course * 1.5)
        elseif r < env.p_init_I_R + env.p_init_I_S
            bstate, variant, abx = I_S, :S, env.abx_course
        elseif r < env.p_init_I_R + env.p_init_I_S + env.p_init_C_R
            bstate, variant, abx = C_R, :R, 0
        elseif r < env.p_init_I_R + env.p_init_I_S + env.p_init_C_R + env.p_init_C_S
            bstate, variant, abx = C_S, :S, 0
        else
            bstate, variant, abx = S, :S, 0
        end

        new_baby = Agent(new_id, BABY, bstate, CLEAN, variant, 0, abx)
        push!(env.babies, new_baby)
    end

    env.tick += 1
end



#---------------------------------------------------------------------------------------------------------------------
# 8.Simulation loop
# Runs the full simulation for a given number of days (142 days here by default)
# and record the number of babies in each state (S, C_S, C_R, I_S, I_R, REC)
function run!(env::Environment; days=142, rng=Random.GLOBAL_RNG)
    record_stats!(env)
    for _ in 1:days
        daily_step!(env, rng)
        record_stats!(env)
    end
end


#---------------------------------------------------------------------------------------------------------------------
# 9. Results and visualisztions
function main()
    println("Starting ESBL-producing Enterobacterales Dynamics in Tanzanian NICU")

    env = init_env()

    days_total = 142
    recoveries_cum = Int[]
    push!(recoveries_cum, 0)

    # record initial stats
    record_stats!(env)

    for _ in 1:days_total
        daily_step!(env, Random.GLOBAL_RNG)
        record_stats!(env)
        push!(recoveries_cum, env.total_recoveries)
    end

    # Compute total infections per day
    infections = []
    for t in 1:length(env.stats[I_S])
        total_today = env.stats[I_S][t] + env.stats[I_R][t]
        push!(infections, total_today)
    end
    peak_inf = maximum(infections)

    # Plotting
    days = 0:days_total

    # Beautify a bit :)
    state_styles = Dict(
        S => ("#5f5f5f", :solid),
        C_S => ("#7262ac", :dot),
        C_R => ("#2e7ebb", :dash),
        I_S => ("#e25508", :dashdot),
        I_R => ("#d92523", :dashdotdot),
        REC => ("#2e974e", :solid)
    )

    plt = plot(
        title="ESBL-producing Enterobacterales Dynamics in NICU (with Total Recoveries)",
        xlabel="Days",
        ylabel="Number of Neonates",
        legend=:right,
        dpi=300,
        size=(1000, 550),
        xtickfont=font(10),
        ytickfont=font(10),
        guidefont=font(13),
        titlefont=font(14),
        legendfontsize=10,
        grid=:y
    )

    # add all states to the plot
    for s in instances(BabyState)
        if s != REC                                 # Exclude REC from the plot
            color, linetype = state_styles[s]
            plot!(plt, days, env.stats[s], label=string(s), color=color, linestyle=linetype, linewidth=2)
        end
    end

    # add Total Recoveries to the plot
    plot!(plt, days, recoveries_cum, label="Total Recoveries", color=:black, linestyle=:dashdot, linewidth=2)

    # make the y-axis all integers
    ymax = maximum(vcat([maximum(env.stats[s]) for s in instances(BabyState)]..., maximum(recoveries_cum)))
    plot!(plt, yticks=0:1:ymax)

    display(plt)

    # Summary output (Some important statistics)
    println("\nPeak infections (I_S + I_R): ", peak_inf)
    println("Total antibiotic days: ", env.total_abx_days)
    println("C_S → C_R selection events:  ", env.selection_events)
    println("Total recoveries (accumulated): ", env.total_recoveries)

    # Average ABX days per baby (based on total simulated population = 165)
   avg_abx_days_per_baby = round(env.total_abx_days / 165, digits=2)
   println("Avg total ABX treatment days per baby: ", avg_abx_days_per_baby)
    
    # Count cumulative resistant infections (I_R)
    total_I_R = sum(env.stats[I_R])
    println("Cumulative resistant infections (I_R): ", total_I_R)

end

main()




#---------------------------------------------------------------------------------------------------------------------
# 10. Run multiple simulations
# Simulate multiple runs of the model and plot the results
function run_multiple!(N_runs::Int=100, days::Int=142)                                      # Define a function to run multiple simulations : times & days
    
    # to obtain the mean and standard deviation of the results
    peak_inf_list = Float64[]
    abx_days_list = Float64[]
    selection_list = Float64[]
    abx_days_per_baby_list = Float64[]

    println("Running $N_runs simulations...")                                              # Annotate a bit

    all_stats = Dict(s => zeros(Float64, days + 1, N_runs) for s in instances(BabyState))  # Create a dictionary to store the stats of each state
    all_infections = zeros(Float64, days + 1, N_runs)                                      # Create a matrix to store the number of infections (I_S + I_R)                                                                                            # for each run
    all_recoveries = zeros(Float64, days + 1, N_runs)                                      # store the cumulative recoveries for each run                                                                                       

    # Run the simulation N_runs times
    for run in 1:N_runs
        env = init_env()

        # run cumulative recoveries
        recoveries_cum = Int[]
        push!(recoveries_cum, 0)

        record_stats!(env)

        for _ in 1:days
            daily_step!(env, Random.GLOBAL_RNG)
            record_stats!(env)
            push!(recoveries_cum, env.total_recoveries)
        end

        for (_, s) in enumerate(instances(BabyState))
            all_stats[s][:, run] = env.stats[s]
        end
        all_infections[:, run] = [env.stats[I_S][t] + env.stats[I_R][t] for t in 1:(days+1)]
        all_recoveries[:, run] = recoveries_cum

        push!(peak_inf_list, maximum(env.stats[I_S] .+ env.stats[I_R]))
        push!(abx_days_list, env.total_abx_days)
        push!(selection_list, env.selection_events)
        push!(abx_days_per_baby_list, round(env.total_abx_days / length(env.babies), digits=2))

    end

    # Create a plot for each state and add the mean & SD
    println("Plotting state trajectories with ± SD...")

    # Beautify a bit:)
    state_styles = Dict(
        S => ("#5f5f5f", :solid),
        C_S => ("#7262ac", :dot),
        C_R => ("#2e7ebb", :dash),
        I_S => ("#e25508", :dashdot),
        I_R => ("#d92523", :dashdotdot),
        REC => ("#2e974e", :solid)
    )

    day_range = 0:days
    _ = ceil(Int, maximum(all_stats[REC])) + 1

    plt1 = plot(title="State-wise Average Dynamics ± SD (Times run=$N_runs)",
        xlabel="Days", ylabel="Number of Neonates", legend=:right, dpi=300,
        size=(800, 400), grid=:y, left_margin=10mm, bottom_margin=8mm, top_margin=6mm,
        guidefont=font(12), xtickfont=font(10), ytickfont=font(10))

    # Use mean and standard deviation to plot the average dynamics of each state
    # Add the ribbon to show the fluctuations
    for s in instances(BabyState)
        μ_run = mean(all_stats[s], dims=2)[:]
        σ_run = std(all_stats[s], dims=2)[:]
        color, style = state_styles[s]
        plot!(plt1, day_range, μ_run, ribbon=σ_run, label=string(s), color=color, linestyle=style, linewidth=2)
    end

    display(plt1)

    # Visualize the total infections (I_S + I_R) over time
    # and add the mean and standard deviation
    println("Plotting total infection trends ± SD...")

    μ_inf = mean(all_infections, dims=2)[:]
    σ_inf = std(all_infections, dims=2)[:]

    plt2 = plot(day_range, μ_inf, ribbon=σ_inf, color=:darkred, lw=2,
        title="Total Number of Infections (sensitive & resistant) ± SD (Times run=$N_runs)",
        xlabel="Days", ylabel="Number of Neonates",
        legend=false, dpi=300, size=(900, 400), grid=:y, left_margin=10mm, bottom_margin=8mm, top_margin=6mm,
        guidefont=font(12), xtickfont=font(10), ytickfont=font(10))

    display(plt2)

    # Visualize the total recoveries (REC) over time
    println("Plotting total recoveries ± SD...")

    μ_rec = mean(all_recoveries, dims=2)[:]
    σ_rec = std(all_recoveries, dims=2)[:]

    plt3 = plot(day_range, μ_rec, ribbon=σ_rec, color=:black, lw=2,
        title="Total Recoveries (accumulated) ± SD (Times run=$N_runs)",
        xlabel="Days", ylabel="Cumulative Recoveries",
        legend=false, dpi=300, size=(900, 400), grid=:y, left_margin=10mm, bottom_margin=8mm, top_margin=6mm,
        guidefont=font(12), xtickfont=font(10), ytickfont=font(10))

    display(plt3)

    # Print summary statistics
    println("\n--- Summary across $N_runs simulations ---")

    # 1. peak infections (I_S + I_R) (average)
    peak_inf_mean = round(mean(peak_inf_list), digits=2)
    peak_inf_sd = round(std(peak_inf_list), digits=2)
    println("Peak infections (mean ± SD): $peak_inf_mean ± $peak_inf_sd")

    # 2. Total ABX days (average)
    abx_days_mean = round(mean(abx_days_list), digits=2)
    abx_days_sd = round(std(abx_days_list), digits=2)
    println("Total ABX days (mean ± SD): $abx_days_mean ± $abx_days_sd")

    # 3. Events of count of C_S → C_R (average)
    selection_mean = round(mean(selection_list), digits=2)
    selection_sd = round(std(selection_list), digits=2)
    println("Selection events (C_S → C_R) (mean ± SD): $selection_mean ± $selection_sd")

    # 4. ABX treatment days per baby (average)
    N_total_babies = 165
    abx_per_baby_mean = round(mean(abx_days_list) / N_total_babies, digits=2)
    abx_per_baby_sd = round(std(abx_days_list) / N_total_babies, digits=2)
    println("Avg total ABX treatment days per baby (mean ± SD): $abx_per_baby_mean ± $abx_per_baby_sd")

    # 5. Cumulative recoveries (at final day)
    final_recoveries = all_recoveries[end, :]  
    rec_mean = round(mean(final_recoveries), digits=2)
    rec_sd = round(std(final_recoveries), digits=2)
    println("Cumulative recoveries at day $days (mean ± SD): $rec_mean ± $rec_sd")

    # 6. Cumulative resistant infections (I_R) across all days
    resistant_inf_per_run = [sum(all_stats[I_R][:, run]) for run in 1:N_runs]
    resistant_inf_mean = round(mean(resistant_inf_per_run), digits=2)
    resistant_inf_sd = round(std(resistant_inf_per_run), digits=2)
    println("Cumulative resistant infections (I_R) (mean ± SD): $resistant_inf_mean ± $resistant_inf_sd")


end

run_multiple!(100)






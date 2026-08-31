using DataFrames
using StatsBase
using Distributions
using Random
using NeutralLandscapes 
using CSV

# Land proportions calculated from Burlington raster data
land_proportions =  [0.2585, 0.2337, 0.1915, 0.1266, 0.0899, 0.0619, 0.0267, 0.0079, 0, 0]

# Data frame of habitat types & movement coefficients (see McClure et al. 2022)
hab_names = ["Deciduous", "DevLo", "Pasture", "DevHi", "Wetlands", "Conifers", "Crops", "Shrub",
                "barrier", "buffer"]
hab_coefs = [0.124, 0, -0.044, -0.496, 0.56, -0.143, -0.556, 0.441, 0, -0.5]

# Burlington:
hab_frame = DataFrame(type = hab_names, prop = land_proportions, coef = hab_coefs)

# Load in functions
include("Functions_sensitivity.jl")

# Create table of param combinations
maxK = [20,25,30]
a_mort = [0.005, 0.0075, 0.01]
j_mort = [0.015, 0.02, 0.025]

all_combos = DataFrame(Iterators.product(maxK, a_mort, j_mort))

#lambda1 = [0.015, 0.02, 0.025, 0.03, 0.035]
#lambda2 = [0.002, 0.003, 0.004, 0.005, 0.006]

#all_combos = DataFrame(Iterators.product(lambda1, lambda2))

# Assign job 
job = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

params = [all_combos[job,1], all_combos[job,2], all_combos[job,3]]

# Simulation function
function the_mega_loop(;years, time_steps, rep, outputs, land_size, maxK, l1, l2, start_cases, amort, jmort)
    # Define average population-level immunity
    seroprev = 0.0

    # create landscape
    landscape = initialize_land(land_size=land_size, barrier_strength = 0)

    # Create array of vaccination probabilities
    vaxprob = fill(seroprev, (land_size,land_size))

    # Fill in buffer zone
    vaxprob[1:5,:] .= 0.6
    vaxprob[(land_size-5):land_size,:] .= 0.6
    vaxprob[:,1:5] .= 0.6
    vaxprob[:,(land_size-5):land_size] .= 0.6

    # Populate landscape
    lil_guys = populate_landscape(seros=seroprev, land_size=land_size)

    # define home coordinates for distance-decay function
    home_coords = deepcopy(lil_guys[:,[1,2,3]])

    for year in 1:years
        for step in 1:time_steps
            # Lots of death
            dont_fear_the_reaper(dat=lil_guys, home=home_coords, k=maxK, amort=amort, jmort=jmort)

            # Move around
            moves = look_around.(lil_guys.x, lil_guys.y, land_size)
            move(moves, lil_guys, home_coords, landscape, 500, -0.001)

            # Spread disease
            spread_disease(dat=lil_guys, home=home_coords, lambda1=l1, lambda2=l2)

            # Immigration as propagule rain
            immigration(dat=lil_guys,home=home_coords,land_size=land_size, disease_rate = 0,
                            sero_rate=0, immigration_rate=1, year=year)  

            # Dispersal
            if 41 <= step <= 45
                # all juveniles go through the dispersal function, but a dispersal distance of 0 is possible
                juvies_leave(lil_guys, home_coords, land_size)
            
                # Not all adults affected by this function, and some have a dispersal distance of 0
                adults_move(lil_guys, home_coords, land_size, year)
            end

            # Reproduction occurs at specific time steps
            if step == 18
                reproduce(lil_guys, home_coords)
            end

            # Vaccine baits are distributed at specific time steps
            if step == 35 
                ORV(dat=lil_guys, land=vaxprob, sero_prob=seroprev)
            end
        
            # Function where some infected guys become symptomatic or recover
            change_state(lil_guys)

            # all guys age 1 week
            lil_guys.age = lil_guys.age .+ 1

            # Update time since infection & disease
            lil_guys.time_since_inf[lil_guys.incubation .== 1] = lil_guys.time_since_inf[lil_guys.incubation.==1] .+ 1
            lil_guys.time_since_disease[lil_guys.infectious .== 1] = lil_guys.time_since_disease[lil_guys.infectious.==1] .+ 1

            # Initialize disease when population stabilizes
            if year == 2 && step == 1
                initialize_disease(dat=lil_guys, nstart=start_cases)
            end       

            elimination = ifelse(sum(lil_guys.incubation) .== 0 .&& sum(lil_guys.infectious) .== 0, "True", "False")

            # Filter out buffer zone
            buffer = filter([:x, :y] => (x, y) -> 5 < x < 55 && 5 < y < 55, lil_guys)
            
            # Calculate summary statistics and append to data frame
            row = [rep, year, step, land_size, maxK, l1, l2, start_cases, amort, jmort, size(buffer,1), 
                    sum(buffer.incubation), sum(buffer.infectious), sum(buffer.vaccinated)/size(buffer,1), 
                    elimination]

            push!(outputs, row)
        end
    end
end

# Run it!
# Create empty data frame

outputs = DataFrame([[], [], [], [], [], [],[],[],[],[],[],[],[],[],[]],
                    ["rep", "year", "week", "land_size", "maxK", "lambda1", "lambda2", "starting_cases",
                    "a_mort", "j_mort", "total_pop", "n_infected", "n_symptomatic", "actual_sero", "elim"])


reps = 20
lam1 = [0.001, 0.00188, 0.00355, 0.00669, 0.0126, 0.0238, 0.0448, 0.0845, 0.159, 0.3]

for rep in 1:(reps-1)
    for i in 1:length(lam1)
        the_mega_loop(years=11, time_steps = 52, rep=rep, outputs = outputs, land_size=60, maxK=25, l1=lam1[i],
                        l2=0, start_cases=0, amort = 0.005, jmort=0.025)

        # Create filename
        filename = string("l1", string(lam1[i]), "rep", string(rep), ".csv")

        # Save results
        CSV.write(filename, outputs)
    end
    println(rep)
end



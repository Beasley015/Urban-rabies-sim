using DataFrames
using StatsBase
using Distributions
using Random
using NeutralLandscapes 
using CSV
#using Dates
#using PProf

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
    vaxprob[56:60,:] .= 0.6
    vaxprob[:,1:5] .= 0.6
    vaxprob[:,56:60] .= 0.6

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
            #spread_disease(dat=lil_guys, home=home_coords, lambda1=l1, lambda2=l2)

            # Immigration can be a propagule rain (steady rate) or a wave (seasonal bursts of high immigration)
            # Remove when doing parameter sensitivity
            #=
            if immigration_type == "propagule"
                immigration(dat=lil_guys,home=home_coords,land_size=land_size, disease_rate = immigration_disease,
                                sero_rate=0, immigration_rate=immigration_rate, year=year)
            elseif immigration_type == "wave"
                if 40 < step < 50
                    immigration(dat=lil_guys,home=home_coords,land_size=land_size, disease_rate = immigration_disease,
                                type="wave", sero_rate=0, immigration_rate=immigration_rate, year=year)
                end
            end
            =#

            # Dispersal
            if step == 43
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
            #=
            if year == 2 && step == 1
                initialize_disease(dat=lil_guys, nstart=start_cases)
            end           
            =#

            elimination = ifelse(sum(lil_guys.incubation) .== 0 .&& sum(lil_guys.infectious) .== 0, "True", "False")

            # Filter out buffer zone
            buffer = filter([:x, :y] => (x, y) -> 5 < x < 55 && 5 < y < 55, lil_guys)
            
            # Calculate summary statistics and append to data frame
            row = [rep, year, step, land_size, maxK, l1, l2, start_cases, amort, jmort, size(buffer,1), 
                    sum(buffer.incubation), sum(buffer.infectious), sum(buffer.vaccinated)/size(buffer,1), 
                    elimination]

            push!(outputs, row)
        end
        println(string("year = ", year))
    end
end

# Run it!
# Create empty data frame

outputs = DataFrame([[], [], [], [], [], [],[],[],[],[],[],[],[],[],[]],
                    ["rep", "year", "week", "land_size", "maxK", "lambda1", "lambda2", "starting_cases",
                    "a_mort", "j_mort", "total_pop", "n_infected", "n_symptomatic", "actual_sero", "elim"])


reps = 10

a_morts = [0.005, 0.0075, 0.01]
j_morts = [0.015, 0.02, 0.025]
Ks = [12, 15, 18]

for i in 1:length(a_morts)
    for j in 1:length(j_morts)
        for k in 1:length(Ks)
            for rep in 1:reps
                the_mega_loop(years=11, time_steps = 52, rep=rep, outputs = outputs, land_size=60, maxK=Ks[k], l1=0.035,
                             l2=0.02, start_cases=10, amort = a_morts[i], jmort=j_morts[j])

        println(string("K = ", Ks[k], ", AdultMortality = ", a_morts[i], ", JuvieMortality = ", j_morts[j], 
                ", rep = ", rep))
            end
        end
    end
end

# Create filename
filename = string("K_mortality_sensitivity.csv")

# Save results
CSV.write(filename, outputs)

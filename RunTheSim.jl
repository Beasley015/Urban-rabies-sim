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
hab_frame = DataFrame(type = hab_names, prop = land_proportions, coef = hab_coefs)

# Load in functions
include("Functions.jl")

# Load in parameters
job = 5#parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
Params = CSV.read("params.csv", DataFrame, skipto=job+1, limit=1, header=1)

# Simulation function
function the_mega_loop(;years, time_steps, seros, rep, immigration_type, immigration_disease, immigration_rate, 
                        outputs)
    # Define average population-level immunity
    seroprev = seros

    # create landscape
    land_size = 60
    landscape = initialize_land(land_size=land_size, barrier_strength = 0, habitats = hab_frame)

    # Create array of vaccination probabilities
    vaxprob = fill(seroprev, (land_size,land_size))

    # Fill in buffer zone
    vaxprob[1:5,:] .= 0.6
    vaxprob[56:60,:] .= 0.6
    vaxprob[:,1:5] .= 0.6
    vaxprob[:,56:60] .= 0.6

    # Populate landscape
    lil_guys = populate_landscape(seros=seroprev)

    # define home coordinates for distance-decay function
    home_coords = deepcopy(lil_guys[:,[1,2,3]])

    for year in 1:years
        for step in 1:time_steps
            # Initialize disease when population stabilizes
            if year == 2 && step == 1
                initialize_disease(lil_guys)
            end

            # Lots of death
            dont_fear_the_reaper(dat=lil_guys, home=home_coords)

            # Move around
            moves = look_around.(lil_guys.x, lil_guys.y, land_size)
            move(moves, lil_guys, home_coords, landscape, 500)

            # Spread disease
            spread_disease(dat=lil_guys, home=home_coords)

            # Immigration can be a propagule rain (steady rate) or a wave (seasonal bursts of high immigration)
            if immigration_type == "propagule"
                immigration(dat=lil_guys,home=home_coords,land_size=land_size, disease_rate = immigration_disease,
                                sero_rate=0.3, immigration_rate=immigration_rate, year=year)
            elseif immigration_type == "wave"
                if 40 < step < 50
                    immigration(dat=lil_guys,home=home_coords,land_size=land_size, disease_rate = immigration_disease,
                                type="wave", sero_rate=0.3, immigration_rate=immigration_rate, year=year)
                end
            end

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

            elimination = ifelse(sum(lil_guys.incubation) .== 0 .&& sum(lil_guys.infectious) .== 0, "True", "False")

            # Filter out buffer zone
            buffer = filter([:x, :y] => (x, y) -> 5 < x < 55 && 5 < y < 55, lil_guys)
            
            # Calculate summary statistics and append to data frame
            row = [rep, year, step, seros, immigration_disease, immigration_rate, immigration_type, size(buffer,1), 
                    sum(buffer.incubation), sum(buffer.infectious), sum(buffer.vaccinated)/size(buffer,1), 
                    elimination]
            push!(outputs, row)
        end
        println(string("Year =", year))
    end
    println(string("Rep = ", rep))
end

# Run it!
# Create empty data frame
outputs = DataFrame([[], [], [], [], [], [],[],[],[],[],[],[]], 
                    ["rep", "year", "week","sero","disease","rate","type", "total_pop", "n_infected", 
                    "n_symptomatic","actual_sero", "elim"])

reps = 50

for rep in 1:reps
    the_mega_loop(years=11, time_steps = 52, seros=Params[!,1][1], rep=rep, immigration_disease = Params[!,3][1], 
                    immigration_type=Params[!,4][1], immigration_rate = Params[!,2][1], outputs = outputs)
end

# Create filename
filename = string("sero",string(Params[!,1][1]),"im_rate",string(Params[!,2][1]),"im_dis",string(Params[!,3][1]),
                    "im_type",string(Params[!,4][1]),".csv")

# Save results
CSV.write(filename, outputs)

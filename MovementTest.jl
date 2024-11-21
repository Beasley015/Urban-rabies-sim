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
hab_frame = DataFrame(type = hab_names, prop = land_proportions, coef = hab_coefs)

# Load in functions
include("Functions.jl")

# Load in parameters
job= 1#parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
Params = CSV.read("params.csv", DataFrame, skipto=job+1, limit=1, header=1)

# Simulation function
function the_mega_loop(;years, time_steps, seros, rep, immigration_type, immigration_disease, immigration_rate, 
                        outputs)

    # create landscape
    land_size = 20
    landscape = initialize_land(land_size=land_size, barrier_strength = 0, habitats = hab_frame)

    # Populate landscape
    lil_guys = populate_landscape(seros=0, guy_density=0.5)

    # define home coordinates for distance-decay function
    home_coords = deepcopy(lil_guys[:,[1,2,3]])

    for year in 1:years
        for step in 1:time_steps
            # Initialize disease at year 5, when population stabilizes
            if year == 5 && step == 1
                initialize_disease(lil_guys)
            end

            # Lots of death
            dont_fear_the_reaper(dat=lil_guys, home=home_coords)

            # Move around
            moves = look_around.(lil_guys.x, lil_guys.y, land_size)
            move(moves, lil_guys, home_coords, landscape, 500, -0.1)

            # Spread disease
            spread_disease(dat=lil_guys, home=home_coords, lambda1=0.015, lambda2=0.0005)

            # Reproduction occurs at specific time steps
            if step == 18
                reproduce(lil_guys, home_coords)
            end
        
            # Function where some infected guys become symptomatic or recover
            change_state(lil_guys)

            # all guys age 1 week
            lil_guys.age = lil_guys.age .+ 1

            # Update time since infection & disease
            lil_guys.time_since_inf[lil_guys.incubation .== 1] = lil_guys.time_since_inf[lil_guys.incubation.==1] .+ 1
            lil_guys.time_since_disease[lil_guys.infectious .== 1] = lil_guys.time_since_disease[lil_guys.infectious.==1] .+ 1

            elimination = ifelse(sum(lil_guys.incubation) .== 0 .&& sum(lil_guys.infectious) .== 0, "True", "False")
            
            if year >= 5
                # get locations of symptomatic guys
                infec = filter(:incubation => x -> x .== 1, lil_guys)

                # put it in a data frame
                frame = DataFrame(year = year, week = step, id = infec.id, x = infec.x, y = infec.y,
                                        inc = infec.incubation, inf = infec.infectious)

                # Calculate summary statistics and append to data frame
                append!(outputs, frame, promote = true)
            end
    
        end
        println(year)
    end
end
# Run it!
# Create empty data frame
outputs = DataFrame([[], [], [], [], [], [], []], 
                    ["year", "week", "id", "x","y", "inc", "inf"])
reps = 1

for rep in 1:reps
    the_mega_loop(years=8, time_steps = 52, seros=Params[!,1][1], rep=rep, immigration_disease = Params[!,3][1], 
                        immigration_type=Params[!,4][1], immigration_rate = Params[!,2][1], outputs = outputs)
end

# Create filename
filename = "mvt_disease.csv"#"c:/users/beasl/documents/urban-rabies-sim/FunctionalityTest/mvt.csv"#string("sero",string(Params[!,1][1]),"im_rate",string(Params[!,2][1]),"im_dis",string(Params[!,3][1]),
                                        #"im_type",string(Params[!,4][1]),".csv")

# Save results
CSV.write(filename, outputs)
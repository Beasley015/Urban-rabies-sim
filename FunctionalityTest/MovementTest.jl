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
include("Functions_smol.jl")

# Load in parameters
job= 1#parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
Params = CSV.read("params.csv", DataFrame, skipto=job+1, limit=1, header=1)

# Simulation function
function the_mega_loop(;years, time_steps, seros, rep, immigration_type, immigration_disease, immigration_rate, 
                        outputs)

    # create landscape
    land_size = 20
    landscape = initialize_land(land_size=land_size, barrier_strength = 0, habitats = hab_frame, movement_test=true)

    # Populate landscape
    lil_guys = populate_landscape(seros=0, guy_density=0.5, movement_test=true)

    # define home coordinates for distance-decay function
    home_coords = deepcopy(lil_guys[:,[1,2,3]])

    for year in 1:years
        for step in 1:time_steps
            # Initialize disease at year 2, when population stabilizes
            if year == 2 && step == 1
                initialize_disease(lil_guys)
            end

            # Lots of death
            dont_fear_the_reaper(dat=lil_guys, home=home_coords, step=step)

            # Move around
            moves = look_around.(lil_guys.x, lil_guys.y, land_size)
            move(moves, lil_guys, home_coords, landscape, 500, -0.001)

            # Spread disease
            spread_disease(dat=lil_guys, home=home_coords)

            # Reproduction occurs at specific time steps
            if step == 18
                reproduce(lil_guys, home_coords)
            end
        
            # Dispersal
            if step == 43
                # all juveniles go through the dispersal function, but a dispersal distance of 0 is possible
                juvies_leave(lil_guys, home_coords, land_size)
                        
                # Not all adults affected by this function, and some have a dispersal distance of 0
                adults_move(lil_guys, home_coords, land_size, year)
            end

            # Function where some infected guys become symptomatic or recover
            change_state(lil_guys)

            # all guys age 1 week
            lil_guys.age = lil_guys.age .+ 1

            # Update time since infection & disease
            lil_guys.time_since_inf[lil_guys.incubation .== 1] = lil_guys.time_since_inf[lil_guys.incubation.==1] .+ 1
            lil_guys.time_since_disease[lil_guys.infectious .== 1] = lil_guys.time_since_disease[lil_guys.infectious.==1] .+ 1

            #elimination = ifelse(sum(lil_guys.incubation) .== 0 .&& sum(lil_guys.infectious) .== 0, "True", "False")

            # Code for testing movement:
            #=
            if year > 1 && step < 43
                df_step = DataFrame(rep=rep, year=year, week=step, id=lil_guys.id, x=lil_guys.x, y=lil_guys.y, 
                                        hab=landscape[CartesianIndex.(lil_guys.x, lil_guys.y)])
                append!(outputs, df_step, promote = true)
            end
            =#
            
            # Code for testing disease transmission:
            #=
            if year >= 2
                # get locations of symptomatic guys
                infec = filter(:incubation => x -> x .== 1, lil_guys)

                # put it in a data frame
                frame = DataFrame(year = year, week = step, id = infec.id, x = infec.x, y = infec.y,
                                        inc = infec.incubation, inf = infec.infectious)

                # Calculate summary statistics and append to data frame
                append!(outputs, frame, promote = true)
            end
            =#
        end
        #println(year)
    end
    # Include this line to save the landscape:
    #CSV.write("ExampleLand.csv",  Tables.table(landscape), writeheader=false)
    
    #println(string("Rep = ", rep))
end

# Run it!
# Create empty data frame
#=
outputs = DataFrame([[], [], [], [], [], [], []], 
                    ["year", "week", "id", "x", "y", "inc", "inf"])
                    =#

# use this for mortality tests
dead_bois = DataFrame([[], [], [], [], [], [],[]], 
            ["step", "n_random_mort", "n_dis_mort", "orphan_mort", "juvie_cc_mort", "adult_cc_mort", "pop_size"])

reps = 1

for rep in 1:reps
    the_mega_loop(years=2, time_steps = 52, seros=Params[!,1][1], rep=rep, immigration_disease = Params[!,3][1], 
                        immigration_type=Params[!,4][1], immigration_rate = Params[!,2][1], outputs = dead_bois)
end

# Create filename
filename = "mvt_mortality_test.csv"
# Save results
CSV.write(filename, dead_bois)

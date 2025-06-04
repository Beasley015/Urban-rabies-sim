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

# Define parameters
Params = DataFrame(seros = 0.0, imm_rate = 0, imm_dis = 0, imm_type = "wave")

# Simulation function
function the_mega_loop(;years, time_steps, seros, rep, immigration_type, immigration_disease, immigration_rate, 
                        outputs, movement_test=false, mortality_test=false)

    # create landscape
    if movement_test==true
        land_size=20
    else
        land_size=60
    end
    
    landscape = initialize_land(land_size=land_size, habitats = hab_frame, movement_test=movement_test)

    # Populate landscape
    lil_guys = populate_landscape(seros=0, guy_density=0.5, movement_test=movement_test)

    # define home coordinates for distance-decay function
    home_coords = deepcopy(lil_guys[:,[1,2,3]])

    for year in 1:years
        for step in 1:time_steps
            # Initialize disease at year 2, when population stabilizes
            if movement_test==false
                if year == 2 && step == 1
                    initialize_disease(lil_guys)
                end
            end

            # Lots of death
            dont_fear_the_reaper(dat=lil_guys, home=home_coords, step=step, mortality_test=mortality_test)

            # Move around
            moves = look_around.(lil_guys.x, lil_guys.y, land_size)
            move(moves, lil_guys, home_coords, landscape, 500, -0.001)

            if movement_test==false
                # Spread disease
                spread_disease(dat=lil_guys, home=home_coords)
            end

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

            if movement_test==false
                elimination = ifelse(sum(lil_guys.incubation) .== 0 .&& sum(lil_guys.infectious) .== 0, "True", "False")
            end

            # Code for testing movement:
            if movement_test==true
                if year > 1 && step < 43
                    df_step = DataFrame(rep=rep, year=year, week=step, id=lil_guys.id, x=lil_guys.x, y=lil_guys.y)
                    append!(outputs, df_step, promote = true)
                end
            end
            
            # Code for testing disease transmission:
            if movement_test==false
                if year >= 2
                    # get locations of symptomatic guys
                    infec = filter(:incubation => x -> x .== 1, lil_guys)

                    # put it in a data frame
                    frame = DataFrame(year = year, week = step, id = infec.id, x = infec.x, y = infec.y,
                                            inc = infec.incubation, inf = infec.infectious)

                    # Calculate summary statistics and append to data frame
                    append!(outputs, frame, promote = true)
                end
            end

            # Code for testing mortality
            if mortality_test==true
                # Code for adding rows to data frame
            end
        end
    end
    # Include this line to save the landscape:
    #CSV.write("ExampleLand.csv",  Tables.table(landscape), writeheader=false)
end

# Run it!
#=
# Data frame for movement tests
outputs = DataFrame([[],[],[],[],[],[]],
                    ["rep", "year", "week", "id", "x", "y"])

# Run it
reps = 5
for rep in 1:reps
    the_mega_loop(years=2, time_steps = 52, seros=Params[!,1][1], rep=rep, immigration_disease = Params[!,3][1], 
                    immigration_type=Params[!,4][1], immigration_rate = Params[!,2][1], outputs = outputs,
                    movement_test=true)
end

# Create filename
filename = "mvt_test.csv"
=#

# Data frame for disease transmission tests
outputs = DataFrame([[], [], [], [], [], [], []], 
                    ["year", "week", "id", "x", "y", "inc", "inf"])

                    

# Data frame for mortality tests
#=
outputs = DataFrame([[], [], [], [], [], [],[]], 
            ["step", "n_random_mort", "n_dis_mort", "orphan_mort", "juvie_cc_mort", "adult_cc_mort", "pop_size"])

reps = 1

for rep in 1:reps
    the_mega_loop(years=2, time_steps = 52, seros=Params[!,1][1], rep=rep, immigration_disease = Params[!,3][1], 
                        immigration_type=Params[!,4][1], immigration_rate = Params[!,2][1], outputs = outputs)
end

# Create filename
filename = "mvt_mortality_test.csv"
=#

# Save results
CSV.write(filename, outputs)

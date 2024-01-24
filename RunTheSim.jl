using DataFrames
using StatsBase
using Distributions
using Random
using NeutralLandscapes # There's an incompatability compiling this with another package, but the ones here aren't it
using PkgCite
using GeoData
using Counters
using ProgressMeter
using CSV

# Land proportions calculated from Burlington raster data
land_proportions =  [0.29950784172189776, 0.22446136010674367, 0.2523743902705777, 0.1341207865782969, 0.06295251219458844, 
                        0.003442920576590765, 0.018693265087385436, 0.0023557976245160445, 0.0020911258394032853, 0.0]

# Data frame of habitat types & movement coefficients
hab_names = ["forest", "developed", "pasture", "wetlands", "herbaceous", "cultivated", "barren", "shrub", "barrier"]
hab_coefs = [0.12, -0.5, -0.03, 0.57, 0, -0.5, 0, 0.4, 0]
hab_frame = DataFrame(type = hab_names, prop = land_proportions[2:10], coef = hab_coefs)

# Load in parameters
########################

# Simulation function
function the_mega_loop(;years, seros, rep, immigration_type, immigration_disease, barrier, outputs)
    # Define weeks per year
    time_steps = 52

    # Define seroprevalence
    seroprev = seros

    # create landscape
    land_size = 60
    landscape = initialize_land(land_size=land_size, barrier_strength = barrier, habitats = hab_frame)

    # Calculate proportion of landscape filled by barrier
    ##############################

    # Populate landscape
    lil_guys = populate_landscape(seros=seroprev)

    # define home coordinates for distance-decay function
    home_coords = deepcopy(lil_guys[:,[1,2,3]])

    for year in 1:years
        for step in 1:time_steps
            # Move around
            moves = look_around.(lil_guys.x, lil_guys.y, land_size)
            move(moves, lil_guys, home_coords, landscape, 500, -0.05)

            # Lots of death
            dont_fear_the_reaper(lil_guys, home_coords, 2)

            # Immigration can be a propagule rain, or a wave at a specific time step
            if immigration_type == "propagule"
                immigration(dat=lil_guys,home=home_coords,land_size=land_size, disease_rate = immigration_disease)
            elseif immigration_type == "wave"
                if step == 30
                    immigration(dat=lil_guys,home=home_coords,land_size=land_size,immigration_rate=200, disease_rate = immigration_disease,
                                type="wave")
                end
            end

            # Juveniles reaching independence (default 30 weeks) disperse
            if size(filter(:age => ==(30), lil_guys),1) > 0 # can change to desired dispersal age
                juvies_leave(lil_guys, home_coords, land_size)
            end

            # Reproduction occurs at specific time steps
            if step == 20 
                reproduce(lil_guys, home_coords)
            end

            # Vaccine baits are distributed at specific time steps
            if step == 45 
                ORV(dat=lil_guys, land_size=land_size, sero_prob=seroprev)
            end 
        
            # Function where some infected guys become symptomatic
            begin_symptoms(lil_guys)

            # all guys age 1 week
            lil_guys.age = lil_guys.age .+ 1

            # Update time since infection & disease
            lil_guys.time_since_inf[lil_guys.incubation .== 1] = lil_guys.time_since_inf[lil_guys.incubation.==1] .+ 1
            lil_guys.time_since_disease[lil_guys.infectious.==1] = lil_guys.time_since_disease[lil_guys.infectious.==1] .+ 1

            # Filter out buffer zone
            buffer = filter([:x, :y] => (x, y) -> 5 < x < 55 && 5 < y < 55, lil_guys)

            # Calculate summary statistics and append to data frame
            row = [rep, year, step, seros, immigration_type, immigration_disease, barrier, barrier_prop, size(buffer,1), sum(buffer.incubation), 
                    sum(buffer.infectious)]
            push!(outputs, row)
            
        end
    end
end

# Run it!
# Create empty data frame
outputs = DataFrame([[], [], [], [], [], [],[],[],[],[],[]], 
                    ["rep", "year", "week","sero","type","rate","barrier_val","barrier_prop","total_pop","n_infected","n_symptomatic"])

reps = 2

p = Progress(length(reps))

for i in 1:length(seros) # Need to create a separate frame for each thread, then combine at the end
    for rep in 1:reps
        the_mega_loop(years=10, seros=seros, rep=rep, immigration_disease = 0.1, immigration_type="wave", barrier = 2, 
                        outputs = outputs)
    end
    next!(p)
end

# Create filename

# Save results
#CSV.write("outputs_full.csv", output)
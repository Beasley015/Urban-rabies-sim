using DataFrames
using StatsBase
using Distributions
using Random
using NeutralLandscapes 
using PkgCite
using CSV

# Land proportions calculated from Burlington raster data
land_proportions =  [0.29950784172189776, 0.22446136010674367, 0.2523743902705777, 0.1341207865782969, 0.06295251219458844, 
                        0.003442920576590765, 0.018693265087385436, 0.0023557976245160445, 0.0020911258394032853, 0.0, 0.0]

# Data frame of habitat types & movement coefficients
hab_names = ["forest", "developed", "pasture", "wetlands", "herbaceous", "cultivated", "barren", "shrub", "barrier", "buffer"]
hab_coefs = [0.12, -0.5, -0.03, 0.57, 0, -0.5, 0, 0.4, 0, -0.5]
hab_frame = DataFrame(type = hab_names, prop = land_proportions[2:11], coef = hab_coefs)

# Load in functions
include("Functions.jl")

# Load in parameters
job= parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
Params = CSV.read("params.csv", DataFrame, skipto=job+1, limit=1, header=1)

# Simulation function
function the_mega_loop(;years, seros, rep, immigration_type, immigration_disease, barrier, outputs)
    # Define weeks per year
    time_steps = 52

    # Define seroprevalence
    seroprev = seros

    # create landscape
    land_size = 60
    landscape = initialize_land(land_size=land_size, barrier_strength = barrier, habitats = hab_frame)

    # Populate landscape
    lil_guys = populate_landscape(seros=seroprev)

    # define home coordinates for distance-decay function
    home_coords = deepcopy(lil_guys[:,[1,2,3]])

    for year in 1:years
        for step in 1:time_steps
            # Initialize disease at year 4, when population stabilizes
            if year == 4 && step == 1
                initialize_disease(lil_guys)
            end

            # Move around
            moves = look_around.(lil_guys.x, lil_guys.y, land_size)
            move(moves, lil_guys, home_coords, landscape, 500, -0.05)

            # Lots of death
            dont_fear_the_reaper(lil_guys, home_coords, 2)

            # Immigration can be a propagule rain (steady rate) or a wave (bursts of high immigration)
            if immigration_type == "propagule"
                immigration(dat=lil_guys,home=home_coords,land_size=land_size, disease_rate = immigration_disease,
                                sero_rate=0.6)
            elseif immigration_type == "wave"
                if year in vcat(2:5, 8:10) && 20 < step < 35
                    immigration(dat=lil_guys,home=home_coords,land_size=land_size, disease_rate = immigration_disease,
                                type="wave", sero_rate=0.6)
                end
            end

            # Juveniles reaching independence (default 40 weeks) disperse
            if size(filter(x -> 20<=x<=75, lil_guys.age),1) > 0 # can change to desired dispersal age
                if step == 40  
                    juvies_leave(lil_guys, home_coords, land_size)
                end
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
            row = [rep, year, step, seros, immigration_disease, barrier, size(buffer,1), sum(buffer.incubation), 
                    sum(buffer.infectious), sum(buffer.vaccinated)/size(buffer,1)]
            push!(outputs, row)
            
        end
    end
end

# Run it!
# Create empty data frame
outputs = DataFrame([[], [], [], [], [], [],[],[],[],[],], 
                    ["rep", "year", "week","sero","rate","barrier", "total_pop", "n_infected", "n_symptomatic","actual_sero"])

reps = 50

for rep in 1:reps
    the_mega_loop(years=14, seros=Params[!,1][j], rep=rep, immigration_disease = Params[!,3][j], 
                    immigration_type=Params[!,4][j], barrier = Params[!,2][j], outputs = outputs)
end

# Create filename
filename = string("sero",string(Params[!,1][1]),"bar",string(Params[!,2][1]),"im_dis",string(Params[!,3][1]),
                                        "im_type",string(Params[!,4][1]),".csv")

# Save results
CSV.write(filename, outputs)

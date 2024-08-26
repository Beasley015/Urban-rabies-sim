# Create simulated landscapes
function initialize_land(;land_size = 60, barrier_strength=0, habitats)
    # 1 cell = 0.5 km by 0.5 km
    # 5-cell buffer in outputs

    # Create neutral landscape based on real habitat proportions
    clustered_land = rand(NearestNeighborCluster(0.2), (land_size, land_size))
    landscape = NeutralLandscapes.classify(clustered_land, land_proportions[2:10])

    # Add buffer habitat
    landscape[1:5,:] .= 10
    landscape[56:60,:] .= 10
    landscape[:,1:5] .= 10
    landscape[:,56:60] .= 10

    #=
    # Create barrier if one exists
    if barrier_strength != 0
        # Pick starting and ending edges for barrier
        edges = ["north", "south", "east", "west"]
        barrier1edge = sample(edges, 2, replace=false)
        barrier2edge = sample(edges, 2, replace=false)

        function create_coord(x)
            y=Vector(undef, 2)

            if x == "north"
                y[1] = rand(1:land_size)
                y[2] = land_size
            elseif x == "east"
                y[1] = land_size
                y[2] = rand(1:land_size)
            elseif x == "south"
                y[1] = rand(1:land_size)
                y[2] = 1
            else 
                y[1] = 1
                y[2] = rand(1:land_size)
            end

            return y
        end
        
        # Start & end coords for barrier
        barrier1start = create_coord(barrier1edge[1])
        barrier1end = create_coord(barrier1edge[2])

        barrier2start = create_coord(barrier2edge[1])
        barrier2end = create_coord(barrier2edge[2])
        
        # Create barriers by calculating slope and intercept from the coords
        function add_barriers(startcoords, endcoords)
            # Calculate slope of barrier line
            if (startcoords[1] - endcoords[1]) == 0
                slope = "undefined"
            else
                slope = (startcoords[2]-endcoords[2])/(startcoords[1]-endcoords[1])
            end

            # Intercept of barrier
            if slope != "undefined"
                intercept = startcoords[2] - slope*startcoords[1]
                # Get the cells the barrier falls on
                x = collect(range(start=minimum([startcoords[1], endcoords[1]]), stop=maximum([startcoords[1],endcoords[1]]), step=0.1))
                y = (slope .* x) .+ intercept
                barrier = DataFrame(x = Int.(round.(x)), y = Int.(round.(y)))
                barrier = unique(barrier)
            else
                x = fill(startcoords[1], land_size)
                y = collect(range(start=1, stop=land_size, step=1))

                barrier = DataFrame(x = x, y = y)
            end

            return(barrier)
        end

        # Build the barrier
        barrier1 = add_barriers(barrier1start, barrier1end)
        barrier2 = add_barriers(barrier2start, barrier2end)
        
        landscape[CartesianIndex.(barrier1.x, barrier1.y)] .= 9.0
        landscape[CartesianIndex.(barrier2.x, barrier2.y)] .= 9.0
        
        # Define barrier strength
        barrier_stren = [1,2,3,4]
        barrier_coefs = [-0.25, -0.5, -1, -10]

        barrier_frame = DataFrame(strength=barrier_stren, coef = barrier_coefs)
        barrier_index = findall(x-> x == barrier_strength, barrier_frame.strength)
        coef = barrier_frame.coef[barrier_index]

        habitats[9,:] = ["barrier", 0.0, coef[1]]
    end
    =#

    return landscape
end

# Initialize raccoon populations
function populate_landscape(;land_size = 60, guy_density = 2.75, seros)
    # Define main area of simulation
    xmin = 6; xmax = 55
    ymin = 6; ymax = 55
    land_area = (xmax-xmin)*(ymax-ymin)

    # get number of guys based on landscape size & density
    nguys = rand(Poisson(land_area*guy_density))

    # Create data frame of guys
    lil_guys = DataFrame(id = string.(collect(1:nguys)), x = convert.(Int64,trunc.(xmax*rand(nguys,1)))[:,1],
        y = convert.(Int64,trunc.(ymax*rand(nguys,1)))[:,1], incubation = 0, time_since_inf = 0, infectious = 0, time_since_disease = 0,
        sex = Int.(rand(Bernoulli(0.5), nguys)), mom = NaN, vaccinated = rand(Bernoulli(seros), nguys), age = rand(52:(52*8), nguys))

    # Remove 0's 
    lil_guys.x[lil_guys.x .== 0] .= 1
    lil_guys.y[lil_guys.y .== 0] .= 1

    # Repeat at a low density to populate the buffer
    xpossible = vcat(1:5, 56:60)
    ypossible = vcat(1:5, 56:60)

    buffer_area = 100

    # Buffer density: typically 4 per km^2, so 1 per cell
    nbuffer = rand(Poisson(buffer_area*1))

    buffer_ids = collect((maximum(parse.(Int64,lil_guys.id))+1):(maximum(parse.(Int64,lil_guys.id))+nbuffer))

    lil_guys_buffer = DataFrame(id = string.(buffer_ids), x = sample(xpossible, nbuffer, replace = true), y = sample(ypossible, nbuffer, replace = true), 
        incubation = 0, time_since_inf = 0, infectious = 0, time_since_disease = 0, sex = Int.(rand(Bernoulli(0.5), nbuffer)),
        mom = NaN, vaccinated = 0, age = rand(52:(52*8), nbuffer))

        lil_guys_buffer.vaccinated[lil_guys_buffer.incubation .!= 1] = rand(Bernoulli(0.6), 
                length(lil_guys_buffer.vaccinated[lil_guys_buffer.incubation .!= 1]))

    lil_guys = [lil_guys; lil_guys_buffer]

    return lil_guys
end

function initialize_disease(dat)
    # Get indices of unvaccinated guys
    unvax = findall(dat.vaccinated .== 0)

    # Choose raccoons to infect
    if length(unvax) > 30
        new_diseases = sample(unvax, 30, replace = false)
    else
        new_diseases = unvax
    end

    # Initialize disease
    dat.incubation[new_diseases] .= 1

    # Start infection timer
    dat.time_since_inf = ifelse.(dat.incubation .== 1, 1, dat.time_since_inf)

    return dat
end

# Create function for the guys to look at their surroundings
function look_around(x,y,land_size)
    # Get vector of all moves
    all_moves = [(x-1, y+1), (x, y+1), (x+1, y+1),
                (x-1, y), (x,y), (x+1, y),
                (x-1, y-1), (x, y-1), (x+1, y-1)]

    # Get indices of impossible moves
    good_moves = findall([(0 .< x[1] .<= land_size) .&& (0 .< x[2] .<= land_size) for x in all_moves])

    # Remove impossible moves
    poss_moves = deepcopy(all_moves[good_moves])

    return poss_moves
end

# Movement and disease function
function move(coords, dat, home, landscape, reso=500, rate=-0.5)
    # Where coords = list of tuples representing possible moves,
    # dat = data frame of agents,
    # reso = width/height of grid cell in meters
    # rate = rate of distance-decay. Based on trial and error so raccoons typically stay ~1km from home

    # Create blank array
    habs = Vector(undef, length(coords))

    # Get habitat type to create weights
    for i in 1:length(coords)
        habs[i] = landscape[CartesianIndex.([x[1] for x in coords[i]], [x[2] for x in coords[i]])]
    end

    # Match habitats to McClure coefficients
    hab_prefs = Vector(undef, length(coords))
    for i in 1:length(coords)
        hab_prefs[i] = hab_frame.coef[convert.(Int64, habs[i])]
    end
   
    # Movement weights as a function of distance from initial coords
    distances = Vector(undef, length(coords))
    dist_weights = Vector(undef, length(coords))
    for i in 1:length(coords)
        home_loc = (home.x[i], home.y[i]) 
        distances[i] = ((([x[1] for x in coords[i]].-home_loc[1]).^2 .+ ([x[2] for x in coords[i]].-home_loc[2]).^2)*reso)/100
        dist_weights[i] = exp.(rate .* distances[i])
    end

    # Combine habitat and distance-decay weights
    weights = Vector(undef, length(dist_weights))
    for i in 1:length(dist_weights)
        weights[i] = dist_weights[i] .* hab_prefs[i]
    end

    # Choose new location
    new_location = sample.(coords, Weights.(weights))

    # Make it a df
    new_spots = DataFrame(x = [x[1] for x in new_location], y = [x[2] for x in new_location])

    # Update data frame
    dat.x = deepcopy(new_spots.x)
    dat.y = deepcopy(new_spots.y)

    # kids follow mom
    kids = findall(dat.age .< 30)

    # Sometimes mom died, so need to filter orphans out to murder later
    poss_indices = [findall(dat.id .== dat.mom[kids][i]) for i in 1:length(kids)]
    blanks = findall(isempty.(poss_indices) .== 1)
    deleteat!(kids, blanks)
    
    indices = vcat(poss_indices...)

    # New kid coords
    x_coords = deepcopy(dat.x[indices])
    y_coords = deepcopy(dat.y[indices])

    dat.x[kids] = deepcopy(x_coords)
    dat.y[kids] = deepcopy(y_coords)
end

function spread_disease(;dat)
    # Find all infected guys
    diseased = filter(:infectious => x -> x .== 1, dat)
    diseased_coords = [(diseased.x[i], diseased.y[i]) for i in 1:nrow(diseased)]

    # Infect raccoons that share cell with diseased guy
        
    if length(diseased_coords) != 0
        # Get all guys with shared coords
        direct_exposure = [intersect(findall(.==(diseased_coords[i][1]), dat.x),findall(.==(diseased_coords[i][2]), dat.y)) 
            for i in 1:length(diseased_coords)]
        direct_exposure = sort(unique(vcat(direct_exposure...)))

        # Remove vaccinated individuals
        direct_exposure = direct_exposure[dat.vaccinated[direct_exposure] .== 0]

        # Infect with set probability
        infections = rand(Bernoulli(0.04), length(direct_exposure))
        direct_exposure = direct_exposure[infections .== 1]

        dat.incubation[direct_exposure] .= 1
    end
end

# Transition from incubation period to infectious period
function begin_symptoms(dat)
    # Probability of transition
    prob = rand.(Beta.(dat.time_since_inf[(dat.incubation .== 1) .& (dat.infectious .== 0)] .+ 0.000000000001,5))
    # If the probability is 0 it does not work; so add a miniscule number

    if length(prob) > 0
        dat.infectious[(dat.incubation .== 1) .& (dat.infectious .== 0)] = reduce(vcat,rand.(Bernoulli.(prob), 1))
    end
end

# Reproduction function
function reproduce(dat, home)
    # Get females of reproductive age
    females = filter([:sex, :age] => (x,y) -> x == 1 && y >= 52, dat)
    #females = filter([:sex, :age] => (x, y) -> x == 1 && y >= 75, dat)
    
    # Choose females that actually reproduce
    reproducing = females[randsubseq(1:size(females,1), 0.9),:]

    # Assign number of offspring to each reproducing female
    noffspring = rand(Poisson(4), size(reproducing,1))

    # Create offspring at location of mother
    devil_spawn = DataFrame(x = Int[], y = Int[], mom = String[])
    for i in 1:length(noffspring)
        for j in 1:noffspring[i]
            push!(devil_spawn, (reproducing.x[i], reproducing.y[i], reproducing.id[i]), promote = true)
        end
    end

    # Fill in missing cols
    devil_spawn.incubation .= 0
    devil_spawn.time_since_inf .= 0
    devil_spawn.infectious .= 0
    devil_spawn.time_since_disease .= 0
    devil_spawn.sex .= Int.(rand(Bernoulli(0.5), size(devil_spawn,1)))
    devil_spawn.vaccinated .= 0
    devil_spawn.age .= 0

    id_vector = collect((maximum(parse.(Int64,dat.id))+1):(maximum(parse.(Int64,dat.id))+size(devil_spawn,1)))
    devil_spawn.id = string.(id_vector)

    # Append to main dataset
    append!(dat, devil_spawn, promote = true)
    append!(home, DataFrame(id = devil_spawn.id, x=devil_spawn.x, y=devil_spawn.y), promote = true)
    
end

# Mortality function
function dont_fear_the_reaper(dat, home, time=2)
    # disease mortality set at 2 weeks, which results in 1 infectious time step
    # due to the order in which functions are executed in the loop
    
    # random mortality
    rand_deaths = rand(Bernoulli(0.0005),size(dat,1))
    deleteat!(dat, findall(rand_deaths==1))
    deleteat!(home, findall(rand_deaths==1))

    # disease mortality
    deleteat!(home, findall(x -> x > time, dat.time_since_disease))
    filter!(:time_since_disease => <=(time), dat)

    # old age mortality
    deleteat!(home, findall(x -> x >= (52*8), dat.age))
    if length(findall(x-> x>=(52*8), dat.age)) > 0
        filter!(:age => <(52*8), dat)
    end
    

    # orphan mortality
    kids = findall(dat.age .< 30)

    no_mom = findall(x -> !(x in dat.id), dat.mom[kids])

    deleteat!(home, no_mom)
    filter!(:mom => !in(dat.mom[kids][no_mom]), dat)

    # Density-related mortality
    # get coordinates where there are multiple guys
    new_location = Vector(undef, size(dat,1))
    for i in 1:size(dat,1)
        new_location[i] = (dat.x[i], dat.y[i]) 
    end

    many_guys = collect(keys(filter(kv -> kv.second > 1, countmap(new_location))))
    indices = [findall(==(x), new_location) for x in many_guys]

    # Find cells with max number of guys or greater
    too_many_guys = findall(length.(indices) .> 10) #can adjust this number

    crowded_spots = many_guys[too_many_guys]

    crowded_indices = Vector(undef, length(crowded_spots))
    for i in 1:length(crowded_spots)
        crowded_indices[i] = findall(x -> x.x == crowded_spots[i][1] && x.y == crowded_spots[i][2], eachrow(dat))
    end

    # Split juveniles and adults for differential mortality
    crowded_indices = sort(unique(vcat(crowded_indices...)))
    crowded_adults = intersect(crowded_indices, findall(x -> x > 52, dat.age))
    crowded_juvies = intersect(crowded_indices, findall(x -> x <= 52, dat.age))

    # Decide who dies
    dead_adults = rand(Bernoulli(0.005), length(crowded_adults))
    dead_juvies = rand(Bernoulli(0.01), length(crowded_juvies))

    dead_guys = sort(vcat(crowded_adults[dead_adults .== 1], crowded_juvies[dead_juvies .== 1]))
    
    if length(dead_guys) > 0
        deleteat!(dat, dead_guys)
        deleteat!(home, dead_guys)
    end  
end

# Vaccination function
function ORV(;dat, land_size, land=nothing, sero_prob)
    if land==nothing 
        # Create array of vaccination probabilities
        land = fill(sero_prob, (land_size,land_size))

        # Fill in buffer zone
        land[1:5,:] .= 0.6
        land[56:60,:] .= 0.6
        land[:,1:5] .= 0.6
        land[:,56:60] .= 0.6
    end

    # Get probs of vaccination at guys' locations
    vaxprob = land[CartesianIndex.(dat.x, dat.y)]

    # Juveniles have lower probs
    kids = findall(dat.age .< 52)
    vaxprob[kids] = vaxprob[kids]./2

    # Remove raccoons that have already been vaccinated
    novax = findall(dat.vaccinated .== 0)
    newvax = vaxprob[novax]

    # Generate vaccine outcomes
    vax = Int.(rand.(Bernoulli.(newvax)))

    # Update data
    dat.vaccinated[novax] = vax

    return dat
end

# Juvenile distribution
function juvies_leave(dat, home, land_size)
    # Get Juveniles
    juvies = dat[findall(x -> 20<x<75, dat.age),:]

    # Create break point so it doesn't get stuck
    niter = 0

    while size(juvies,1) > 0
        niter = niter + 1

        # Pick a direction from list of inline functions
        upleft(x,y)=[x-1, y+1]; up(x,y)=[x, y+1]; upright(x,y)=[x+1, y+1]
        left(x,y)=[x-1, y]; right(x,y)=[x+1, y]
        downleft(x,y)=[x-1, y-1]; down(x,y)=[x, y-1]; downright(x,y)=[x+1, y-1]
        directions = rand([upleft, up, upright, left, right, 
                        downleft, down, downright], size(juvies,1))

        # Get dispersal distance
        distances = rand(Poisson(5), size(juvies,1))

        # RUN!
        coords = Vector(undef, size(juvies,1))

        for i in 1:length(distances)
            coords[i] = [juvies.x[i], juvies.y[i]]
                for j in 1:distances[i]
                coords[i] = directions[i](coords[i][1], coords[i][2])
            end
        end

        # Remove juveniles that left the landscape
        new_juvies = juvies[all.(x-> 0<x<land_size, coords),:]
        gone_juvies = juvies[all.(x-> x<0 || x > land_size, coords),:]

        # Update full data frame
        indices=findall(x-> x in(new_juvies.id), dat.id)
        gone_indices = findall(x-> x in(gone_juvies.id), dat.id)
        dat[indices,:]=new_juvies
        deleteat!(dat, gone_indices)

        # Update home coords data frame
        home.x[indices] = new_juvies.x
        home.y[indices] = new_juvies.y
        home.id[indices] = new_juvies.id
        deleteat!(home, gone_indices)
    
        # Find coordinates with multiple guys
        new_location = Vector(undef, size(dat,1))
        for i in 1:size(dat,1)
            new_location[i] = (dat.x[i], dat.y[i]) 
        end
    
        many_guys = collect(keys(filter(kv -> kv.second > 1, countmap(new_location))))
        indices = [findall(==(x), new_location) for x in many_guys]
    
        # Find cells with less than max number of guys
        enough_guys = findall(length.(indices) .<= 8) #can adjust this number
    
        good_spots = many_guys[enough_guys]
    
        good_indices = Vector(undef, length(good_spots))
        for i in 1:length(good_spots)
            good_indices[i] = intersect(findall(x -> x == good_spots[i][1], new_juvies.x), findall(x -> x == good_spots[i][2], new_juvies.y))
        end    

        deleteat!(new_juvies,sort(unique(vcat(good_indices...))))

        # Start the cycle again
        juvies = deepcopy(new_juvies)

        if niter > 5
            break
        end
    end
end

# Immigration function
function immigration(;dat, home, land_size, immigration_rate=5, sero_rate=0, disease_rate=0.3, type="propagule")
    if type == "wave"
        immigration_rate = immigration_rate*20
    end

    # get number of immigrants
    n_new = rand(Poisson(immigration_rate))

    # Data frame of immigrants
    immigrants = DataFrame(id = string.(collect(range(start=maximum(parse.(Int, dat.id)),stop=maximum(parse.(Int, dat.id))+(n_new-1),step=1))), 
                        x = 0, y = 0, incubation = 0, time_since_inf = 0, infectious = 0, time_since_disease = 0, 
                        sex = Int.(rand(Bernoulli(0.5), n_new)), mom = NaN, vaccinated = 0, age = rand(Poisson(65), n_new))

    # Initialize disease
    immigrants.incubation = rand(Bernoulli(disease_rate), n_new)
    immigrants.time_since_inf = ifelse.(immigrants.incubation .== 1, 1, immigrants.time_since_inf)

    # Initialize immunity
    immigrants.vaccinated[immigrants.incubation .!= 1] = rand(Bernoulli(sero_rate), length(immigrants.vaccinated[immigrants.incubation .!= 1]))

    # Get starting edges
    edges = ["north", "east", "south", "west"]

    immigrant_edges = sample(edges, n_new)

    for i in 1:length(immigrant_edges)
        if immigrant_edges[i] == "north"
            immigrants.x[i] = rand(1:land_size)
            immigrants.y[i] = land_size
        elseif immigrant_edges[i] == "east"
            immigrants.x[i] = land_size
            immigrants.y[i] = rand(1:land_size)
        elseif immigrant_edges[i] == "south"
            immigrants.x[i] = rand(1:land_size)
            immigrants.y[i] = 1
        else 
            immigrants.x[i] = 1
            immigrants.y[i] = rand(1:land_size)
        end
    end

    # Pick a direction from list of inline functions
    upleft(x,y)=[x-1, y+1]; up(x,y)=[x, y+1]; upright(x,y)=[x+1, y+1]
    left(x,y)=[x-1, y]; right(x,y)=[x+1, y]
    downleft(x,y)=[x-1, y-1]; down(x,y)=[x, y-1]; downright(x,y)=[x+1, y-1]
    directions = rand([upleft, up, upright, left, right, 
                    downleft, down, downright], size(immigrants,1))

    # Get dispersal distance
    distances = rand(Poisson(10), size(immigrants,1))

    # RUN!
    coords = Vector(undef, size(immigrants,1))

    for i in 1:length(distances)
        coords[i] = [immigrants.x[i], immigrants.y[i]]
        for j in 1:distances[i]
            coords[i] = directions[i](coords[i][1], coords[i][2])
        end
    end

    # Keep immigrants that stayed in the landscape
    immigrants.x = [coords[i][1] for i in 1:length(coords)]
    immigrants.y = [coords[i][2] for i in 1:length(coords)]
    filter!([:x, :y] => (x,y) -> 0 < x <= land_size && 0 < y <= land_size, immigrants)

    # Append immigrants to main dataset & home coords
    dat = [dat;immigrants]
    home = [home;DataFrame(id = immigrants.id, x = immigrants.x, y = immigrants.y)]
end
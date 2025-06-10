#Necessary packages. They must be previously installed for the code to work.
#Note that LightGraphs is outdated, but works. Graphs is more recent and has essentially the same functions.
using SparseArrays, LinearAlgebra, StatsBase, Statistics, LightGraphs, PoissonRandom, DelimitedFiles, Random, Dates, Distributions

#set working directory, where files will be saved
cd("/home")

#Creates a seed for generating random numbers based on current time
seed = Dates.value(convert(Dates.Millisecond, Dates.now()))
Random.seed!(seed)

#This function rewires the network
#Na is the number of nodes, La the side length of the lattice, and pa the rewiring probability
function adjacency(Na::Int64, La::Int64, pa::Float64)

        #Generate a regular square latice with periodic boundary conditions and size La x La
        g::SimpleGraph{Int64} = LightGraphs.grid(Int[La,La], periodic=true)
        #Collect all edges present in the regular lattice
        edge_list::Vector{LightGraphs.SimpleGraphs.SimpleEdge{Int64}} = collect(edges(g))

        #The number of edges to be rewired, drawn from a binomial distribution with 2*Na trials (the total number of edges) and
        #success (rewiring) probability pa.
        number_rewire::Int64 = rand(Binomial(2*Na, pa))
        #Choose which elements from edge_list will be rewired, without replacement since an edge cannot be rewired twice.
        rewire::Vector{Int64} = sample(1:length(edge_list), number_rewire, replace=false) #edges to be rewired

        #Loop to rewire the edges listed in rewire
        for i = 1:length(rewire)
                #This is the node with lowest index in the edge. Nodes are indexed such that the top right is 1, and the index
                #grows the right and down.
                node1::Int64 = src(edge_list[rewire[i]])
                #This is the node in the edge with the highest index
                node2::Int64 = dst(edge_list[rewire[i]])

                #Chooses a new node at random
                new::Int64 = rand(1:Na)
                #Condition to not recreate the edge that will be removed, nor introduce a double or self edge. While any of these
                #happens choose a new random node.
                while new == node1 || new == node2 || has_edge(g, node1, new)
                        new = rand(1:Na)
                end

                #Remove the previous edge
                rem_edge!(g, node1, node2)
                #Replace with the new edge
                #These operations only need to be done once since the graph is undirected and Julia understands the adjacency matrix
                #must remain symmetric
                add_edge!(g, node1, new)
        end

        #Return the rewired graph
        return g
end

#Performs spin flip trials and updates the lattice
function metropolis(Nm, Lm, probm, total_time, spinsm, bm, T, seed)
	
	#Initializes empty arrays that will keep the order parameter and population at each site at each time.
        #They will be dynamically allocated
	order = Float64[]
	spindist = Float64[]
	#Keeps values of the exponential of the reciprocal of temperature, for the purpose of approving spin flips that increase energy
	expo = exp.(bm*(-1:-1:-100))

	#Spin updates for total_time time steps
	for i = 1:total_time
		#Get the rewired network at a given time
		g = adjacency(Nm, Lm, probm)

		#choose which spins to test, at each time step we try to flip N spins (N = total number of nodes), chosen WITH replacement
		test = rand(1:Nm, Nm)
		#Try to flip each one of the spins chosen above
		for j = 1:Nm
			#Get the sum of the neighboring spins
			nn = sum(spinsm[neighbors(g, test[j])])
			#Calculate the energy difference if the spin flips
			delta = 2*spinsm[test[j]]*nn
			#If the flip decreases the energy, accept it
			if delta <= 0 
				spinsm[test[j]] = -1*spinsm[test[j]]
			#If the flip increases the energy, flip following the Gibbs probability
			elseif rand() < expo[delta]
				spinsm[test[j]] = -1*spinsm[test[j]]
			end
		end
		
		#Appends the average spin (aka magnetization) in the order parameter array
		append!(order, mean(spinsm))
		#Appends the spin values in the spins array
		if i%5000 == 0 #saves individual spins every 5000 times steps
                        append!(spindist, spinsm)
                end
	end
	
	#This is to filter which order parameter values we save, since we don't need a very fine time series
	fin_order = Float64[]
        i = 1
        while i < length(order)
                append!(fin_order, order[i])
                if i < 10^7 #finer data at the beginning, times are 1, 2, 3, ..., 10000, 10010, ..., 100000.
                        i += 10
                else
                        i += 500
                end
        end
	#Write files
        writedlm("annealising-order-p$probm-L$Lm-T$T.csv", [fin_order], ',')
        writedlm("annealising-spins-p$probm-L$Lm-T$T.csv", [spindist], ',')
end

#This function selects the parameter values to run the simulation
function sim_run(linear_size, area, probability, bet, total_time, T, seed)

        for i = 1:length(linear_size)
                L = linear_size[i]
                N = area[i]
                for j = 1:length(sig)
                        p = probability
                        b = bet[j]
			tp = T[j]
			#Initial spin distribution is random
                        spins =  sign.(rand(Int, N))
			metropolis(N, L, p, total_time, spins, b, tp, seed)
                end
        end
end

#System sizes
linear_size::Vector{Int64} = [16, 32, 64, 128]
area::Vector{Int64} = linear_size.^2
#Rewiring probabilities
probability::Vector{Float64} = ones(1)
append!(probability, range(0.0, 0.2, length=11), 0.6)
#Temperatures
temp = range(3, 2, length=21)
#Reciprocal of temperature
beta = 1 ./temp
#Run jobs in parallel so each job has a specific temperature
sig = beta[parse(Int, ARGS[1])]
T = temp[parse(Int, ARGS[1])]

#The code commented below was to sample over temperatures closer to the critical temperature.
#nc = [2.75, 2.27, 2.507, 2.59, 2.645, 2.684, 2.713, 2.738, 2.756, 2.77, 2.783, 2.792, 2.787]
#temp = ones(13,2)
#for i=1:length(nc)
#        temp[i,:] = range(nc[i]-0.06, nc[i]+0.06, 2)
#end
#beta = 1 ./temp
#sig = beta[parse(Int, ARGS[1]),:]
#T = temp[parse(Int, ARGS[1]),:]
#probab = probability[parse(Int, ARGS[1])]

global total_time = 5*10^6

#Run the simulation
sim_run(linear_size, area, probab, sig, total_time, T, seed)





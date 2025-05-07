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

#This function updates the population according to the noisy population growth equation and dispersal
#Ni is the number of nodes, Li the side length of the lattice, probi the rewiring probability, total_time the number of time steps
#to be considered, popi the array with the initial population, eps the dispersal fraction, sigmai the noise level, r the growth rate,
#and seed the random number seed.
function iter(Ni::Int64, Li::Int64, probi::Float64, total_time::Int64, popi::Vector{Float64}, eps::Float64, sigmai::Float64, r::Float64, seed::Int64)
	
	#Initializes empty arrays that will keep the order parameter and population at each site at each time. 
	#They will be dynamically allocated
	order = Float64[] 
	populationaux = Float64[]

	#Copy the initial populations
	aux_pop = copy(popi)
	#Population update over total_time time steps
	for k = 1:total_time
		#Get the rewired network at a given time
		g = adjacency(Ni, Li, probi)
		#Calculate the new population following a noisy Ricker growth
		popi2 = zeros(Ni)
		popi2 = @. aux_pop*exp(sigmai*randn()+r*(1-aux_pop))#ricker growth with noise
		#Copies the new popualtion to account for dispersal
		popi3 = copy(popi2)
		#Accounts for dispersal in each of the nodes
		for i = 1:Ni
			#The positive term is for individuals coming into node i from each of its neighbors
			#The negative term is for individuals leaving node i to each of its neighbors
       			popi2[i] = popi3[i] + eps*sum(popi3[neighbors(g, i)]) - eps*popi3[i]*length(neighbors(g,i)) #dispersal
       		end
		
		#Calculate the two-cycle variable by taking the first difference in population
		#popi2 is the current population, aux_pop is the population before dispersal AND noisy local growth
		op = (-1)^k*(popi2-aux_pop)/2

		#Appends the average of the two-cycle variable in the order parameter array
		append!(order, mean(op))
		#Appends the population (before and after local growth+dispersal) in the population array
		if k%5000 == 1 #saves population every 5000 times steps
			append!(populationaux, transpose(aux_pop), transpose(popi2))
		end
		#copies the updated population to be the starting population at the next time step
		aux_pop = copy(popi2)
	end

	#This is to filter which order parameter values we save, since we don't need a very fine time series
	fin_order = Float64[]
	i = 1
	while i < length(order)
		append!(fin_order, order[i])
		if i < 2*10^7 
			i += 10
		else
			i += 500
		end
	end
	
	#Write files
	writedlm("order-p$probi-L$Li-noise$sigmai.csv", [fin_order], ',')
        writedlm("pop-p$probi-L$Li-noise$sigmai.csv", [populationaux], ',')

end

#This function runs the noisy ricker map without dispersal for 100 time steps. 
#It can be used to ensure initial conditions will be random and in the equilibrium state of local dynamics for every node before dispersal is turned on. These do not count as simulation time steps.
function ricker100(x::Vector{Float64}, r::Float64, sig::Float64)
	for i=1:100
		x = @. x*exp(sig*randn()+r*(1-x))::Vector{Float64}
	end
	return x

end

#This function selects the parameter values to run the simulation
function sim_run(linear_size::Vector{Int64}, area::Vector{Int64}, probability::Vector{Float64}, noise::Float64, total_time::Int64, eps::Float64, r::Float64)#, repeat)

	for i = 1:length(linear_size)
        	L::Int64 = linear_size[i]
        	N::Int64 = area[i]
        	for j = 1:length(probability)
			p::Float64 = probability[j]
			sigma::Float64 = noise
			#pop::Vector{Float64} = ricker100((1.72*rand(N).+0.14), 2.3, sigma) #initial pop is node equilibrium with given value of noise, rand(N) modified to give equal probability for each phase of oscillation
			pop = 0.2*randn(N).+1.6 #Initial population is ordered with some standard deviation
			iter(N, L, p, total_time, pop, eps, sigma, r, seed)
        	end
	end
end

#System sizes
linear_size::Vector{Int64} = [16, 32, 64, 128]
area::Vector{Int64} = linear_size.^2
#Rewiring probabilities
probability::Vector{Float64} = ones(1)
append!(probability, range(0.0, 0.2, length=11), 0.6)
#Noise values
noise = range(0.14, 0.195, length=12)
#Running 12 jobs in parallel, each for one noise level
sig = noise[parse(Int, ARGS[1])]

#Total simulation time
const total_time::Int64 = 10^7
#Dispersal fraction
const eps::Float64 = 0.025
#Local growth rate
const r::Float64 = 2.3

#Runs the simulation
sim_run(linear_size, area, probability, sig, total_time, eps, r)#, repeat)

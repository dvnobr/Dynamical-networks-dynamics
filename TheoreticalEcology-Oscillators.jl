import Pkg
#Pkg.add("CSV")
#Pkg.add("DataFrames")
#Pkg.add("Distributions")
#Pkg.add("Tables")
#Pkg.add("PoissonRandom")
#Pkg.add("LightGraphs")
#Pkg.add("DelimitedFiles")
using SparseArrays, LinearAlgebra, StatsBase, Statistics, LightGraphs, PoissonRandom, DelimitedFiles, Random, Dates, Distributions#, DataFrames, CSV, Tables, Dates

cd("/home/dvnobr/Data-p2")

#Random.seed!(seed) - time dependent random seed
seed = Dates.value(convert(Dates.Millisecond, Dates.now()))
Random.seed!(seed)

#This function rewires the network
function adjacency(Na::Int64, La::Int64, pa::Float64)
	#g starts as a regular 2D lattice with periodic boundary conditions
        g::SimpleGraph{Int64} = LightGraphs.grid(Int[La,La], periodic=true)
	
	#collects all edges in g
        edge_list::Vector{LightGraphs.SimpleGraphs.SimpleEdge{Int64}} = collect(edges(g))

        #Chooses how many edges to rewires, from a total of 2*N with probability p
	number_rewire::Int64 = rand(Binomial(2*Na, pa)) #number of edges to be rewired, binomial process with mean p*N
        #Chooses the edges from edge_list that will be rewired
	rewire::Vector{Int64} = sample(1:length(edge_list), number_rewire, replace=false) #edges to be rewired

	#Rewiring the chosen edges
        for i = 1:length(rewire)
		#Collects the two nodes connected by the selected edge
                node1::Int64 = src(edge_list[rewire[i]])
                node2::Int64 = dst(edge_list[rewire[i]])
		
		#Choose a new node
                new::Int64 = rand(1:Na)
		#Choose one of the nodes from the edge with probability 0.5
                nold = rand([node1, node2])
		#Confirms that connecting nold to new is not creating a self-edge, recreating the previous edge, or adding an already existing edge.
                while new == node1 || new == node2 || has_edge(g, nold, new)
                        new = rand(1:Na)
                end
		#Removes the edge between nodes 1 and 2
                rem_edge!(g, node1, node2)
		#Adds the edge between nold and new
                add_edge!(g, nold, new)
        end
	
	#Returns the rewired network
        return g
end

function iter(Ni::Int64, Li::Int64, probi::Float64, total_time::Int64, popi::Vector{Float64}, eps::Float64, sigmai::Float64, r::Float64, seed::Int64, rep)
	
	#Arrays containing the order parameter values over time and the population abundance at each site will be allocated dynamically
        order = Float64[] 
	populationaux = Float64[]

        aux_pop = copy(popi)
        for k = 1:total_time
                #Gets the rewired network
		g = adjacency(Ni, Li, probi) #keeps the adjacency matrix and degree of each node

                popi2 = zeros(Ni)
		#Performs noisy Ricker growth
                popi2 = @. aux_pop*exp(sigmai*randn()+r*(1-aux_pop))#ricker growth with noise
                popi3 = copy(popi2)
		#Dispersal
                for i = 1:Ni
                        popi2[i] = popi3[i] + eps*sum(popi3[neighbors(g, i)]) - eps*popi3[i]*length(neighbors(g,i)) #dispersal
                end
		#Calculates the two-cycle variable at each site
                op = (-1)^k*(popi2-aux_pop)/2
		#Appends the mean two-cycle variable to the order parameter array
                append!(order, mean(op))
		if k%5000 == 1 #in [1,11,101,1001] #saves population every 5000 times steps
                        append!(populationaux, transpose(aux_pop), transpose(popi2)) #Save population at 2 time steps in order to be able to determine first difference
                end
                aux_pop = copy(popi2)
        end
        #saving files
        fin_order = Float64[]
        i = 1
	#This is to save data only at specific time steps (to save memory)
        while i < length(order)
               append!(fin_order, order[i])
               if i < 2*10^7 #finer data at the beginning
                       i += 10
               else
                       i += 500
               end
        end
	#Save files
        writedlm("initial-order-p$probi-L$Li-noise$sigmai.csv", [fin_order], ',')
        writedlm("initial2-pop-p$probi-L$Li-noise$sigmai-rep$rep.csv", [populationaux], ',')

end

function sim_run(linear_size::Vector{Int64}, area::Vector{Int64}, probability::Vector{Float64}, noise::Float64, total_time::Int64, eps::Float64, r::Float64, repeat)
	#SELECTS PARAMETER VALUES
        for i = 1:length(linear_size)
                L::Int64 = linear_size[i]
                N::Int64 = area[i]
                for j = 1:length(probability)
                        for k = 1:repeat
                                p::Float64 = probability[j]
                                sigma::Float64 = noise
                                #pop::Vector{Float64} = ricker100((1.72*rand(N).+0.14), 2.3, sigma) #initial pop is node eq with given value of noise, rand(N) modified to give equal probability for each phase of oscillation
                                #pop = 0.2*randn(N).+1.6
                                pop = 0.2*randn(N).+1.6 #readdlm("paper2-pop-p$p-L$L-noise$sigma.csv",',')[end-N+1:end]
                                iter(N, L, p, total_time, pop, eps, sigma, r, seed, k)
                        end
                end
        end
end
linear_size::Vector{Int64} = [16, 32, 64, 128]
area::Vector{Int64} = linear_size.^2
probability::Vector{Float64} = ones(1)
append!(probability, range(0.0, 0.2, length=11), 0.6)

noise = range(0.14, 0.195, length=12)
#Running in parallel, one job for each noise value
sig = noise[parse(Int, ARGS[1])]

const total_time::Int64 = 2*10^6
const eps::Float64 = 0.025
const r::Float64 = 2.3
const repeat = 1#00 Option to change the number of different realizations. This is only relevant for out-of-equilibrium measurements

sim_run(linear_size, area, probability, sig, total_time, eps, r, repeat)

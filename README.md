# Dynamical networks dynamics
This repository includes the Julia codes for two different systems in dynamical networks: Ising spins and ecological oscillators. Nodes on the network interact with coupled nodes through edges that mostly connect nearest neighbors, but there is a dynamical rewiring on the network to account for some franction of long-range coupling. The methods and results are on my paper - Behavior of Ising spins and ecological oscillators on dynamically rewired small-world networks, which is on arXiv and has been accepted for publication in Physical Review E.

oscillators-PRE.jl is the code for ecological oscillators, which are populations that oscillate according to a noisy Ricker map in the period-2 cyclic regime. Different nodes interact via migration. This code was used for the paper submitted to PRE.

ising-PRE.jl is the code for Ising spins, which can have a value of +1 or -1. Spins tend to align with neighbors in a description of ferromagnetism.

TheoreticalEcology-Oscillators.jl is similar to oscillators-PRE.jl, with very minor changes in the rewiring process and different initial conditions (ordered, instead of random). This code was used in the paper submitted to Theoretical Ecology.

To cite this, please use https://doi.org/10.5281/zenodo.15466486.

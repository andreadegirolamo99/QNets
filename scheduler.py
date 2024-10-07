import os

N = 10
topologies = ["square_lattice", "diagonal_square_lattice", "honeycomb", "fully_connected_honeycomb"]
sigmas = [0, 0.01, 0.07]

for topology in topologies:
    for sigma in sigmas:
        if sigma == 0:
            topology = "noisy_" + topology
        os.system('sbatch ./qnets.sh ' + topology + ' ' + str(N) + ' ' + str(sigma))
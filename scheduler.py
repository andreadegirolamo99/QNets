import os
import numpy as np

N = 10
topologies = ["square_lattice"]#["square_lattice", "diagonal_square_lattice", "honeycomb", "fully_connected_honeycomb"]
sigmas = [0]#[0, 0.01, 0.07]
lambdas = np.arange(0.5, 1.005, 1/200)
n_samples = 1
samples = []

for topology in topologies:
    if not os.path.isdir("out_files/{}_{}".format(topology, N)):
        os.mkdir("out_files/{}_{}".format(topology, N))

    for sigma in sigmas:
        if sigma == 0:
            if not os.path.isdir("out_files/{}_{}/mean".format(topology, N)):
                os.mkdir("out_files/{}_{}/mean".format(topology, N))
        else:
            if not os.path.isdir("out_files/{}_{}/std_dev={}".format(topology, N, sigma)):
                os.mkdir("out_files/{}_{}/std_dev={}".format(topology, N, sigma))

            s_num = 1 if not os.listdir("out_files/{}_{}/std_dev={}".format(topology, N, sigma)) else max([int(sample) for sample in os.listdir("out_files/{}_{}/std_dev={}".format(topology, N, sigma))]) + 1
            for sample in range(n_samples):
                os.mkdir("out_files/{}_{}/std_dev={}/{}".format(topology, N, sigma, s_num))
                samples.append(s_num)
                s_num += 1



for topology in topologies:
    for (i, lambda_) in enumerate(lambdas):
        for sigma in sigmas:
            if sigma > 0:
                topology_ = "noisy_" + topology
                for sample in samples:
                    os.system('./qnets.sh ' + topology_ + ' ' + str(N) + ' ' + str(i+1) + ' ' + str(lambda_) + ' ' + str(sigma) + ' ' + str(sample))
            else:
                os.system('./qnets.sh ' + topology + ' ' + str(N) + ' ' + str(i+1) + ' ' + str(lambda_) + ' ' + str(sigma) + ' ' + str(0))
# TODO:
# - fazer com que os plots a 1D mostrem a label (isto eh algo especifico do getdist, nao sei como alterar)

# imports
from getdist import plots, MCSamples
import matplotlib.pyplot as plt
import numpy as np
import pandas
import os

# specify the output folders with the output of Stan
folders = ["output3/cc", "output3/snia", "output3/bao", "output3/together"]
legends = ["CC", "SnIa", "BAO", "Combined"]

# specify the parameters that were constrained in each model
paramsperfolder = [["H0", "Om"], [ "Om", ],["H0", "Om"], ["H0", "Om"]]
labelsperfolder = [["H0", "\\Omega_m" ], [ "\\Omega_m"], ["H0", "\\Omega_m" ], ["H0", "\\Omega_m" ]]

# get 'MCSamples' object for each run
mcsamples = []
for i in range(0, len(folders)):
    folder  = folders[i]
    params  = paramsperfolder[i]
    labels  = labelsperfolder[i]
    ndim    = len(params)
    legend  = legends[i]

    # get chains for each run
    files = sorted(os.popen(f"find {folder} -type f -name '*.csv'").read().split())
    chains = len(files)

    # get the samples from each chain
    samples = len(pandas.read_csv(files[0], comment="#")[params[0]])
    flatsamples = np.empty([samples*chains, ndim])
    chainN = 0
    for file in files:
        data = pandas.read_csv(file, comment="#")

        for i in range(0, len(params)):
            flatsamples[chainN::chains, i] = data[params[i]]

        chainN += 1

    mcsamples.append(MCSamples(samples=flatsamples, names=params, labels=labels, label=legend))

# make the corner plot
g = plots.get_subplot_plotter()
g.triangle_plot(mcsamples, filled=True)
plt.show()
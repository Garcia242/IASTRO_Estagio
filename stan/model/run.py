# imports
from getdist import plots, MCSamples
from cmdstanpy import CmdStanModel
import matplotlib.pyplot as plt
from numpy.random import normal
import arviz as az

# create a model instance
model = CmdStanModel(
    #compile="force",            # force model recompilation
    stan_file="model/cc.stan",  # Stan model file location
    cpp_options={
        #"STAN_NO_RANGE_CHECKS": "TRUE",  # don't check for elements out of bounds
        "STAN_THREADS": "TRUE",          # run multiple chains in parallel
        "STAN_CPP_OPTIMS": "TRUE"        # optimizations recommended by the Stan development team
    }
)

# configure and fit the model
fit = model.sample(
    #show_console=True,            # useful for debugging
    output_dir="output/cc",       # save output to specified folder
    data="data/cc.json",          # the location of the data file
    iter_warmup=500,              # the number of warmup steps
    iter_sampling=500,            # the number of sampling steps
    inits={"h": 0.7, "Om": 0.3},  # initial values for each parameter (should be randomly chosen for each chain)
    chains=4,                     # number of chains to run in total
    parallel_chains=4             # number of chains to run in parallel (if possible)
)

# print information provided by very useful utilities available in Stan
print("============================== fit summary ==============================")
print(fit.summary())
print("============================ fit diagnostics ============================")
print(fit.diagnose())

# show the traceplot
posterior = az.from_cmdstanpy(posterior=fit)
az.plot_trace(posterior, var_names=("h", "Om"), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show corner plot
samples = samples = [fit.stan_variable("h"), fit.stan_variable("Om")]
mcsamples = MCSamples(samples=samples, names=["h", "Om"], labels=["h", "\\Omega_m"])
g = plots.get_subplot_plotter()
g.triangle_plot(mcsamples, filled=True)
plt.show()

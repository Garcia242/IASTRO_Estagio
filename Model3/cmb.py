# imports
from getdist import plots, MCSamples
from cmdstanpy import CmdStanModel
import matplotlib.pyplot as plt
from numpy.random import normal
import arviz as az

# create a model instance
model = CmdStanModel(
    stan_file="Model3/Model/cmb3_1.stan",                         # Stan model file location
    cpp_options={
        #"STAN_NO_RANGE_CHECKS": "TRUE",  # don't check for elements out of bounds
        "STAN_THREADS": "TRUE",          # run multiple chains in parallel
        "STAN_CPP_OPTIMS": "TRUE",        # optimizations recommended by the Stan development team
        "force_one_process_per_chain": "True"
    }
)
# configure and fit the model
fit = model.sample(
    data="Model3/Data/cmb.json",                                # the location of the data file
    #output_dir="output3/cmb",
    iter_sampling=500,                                                 # the number of sampling steps
    iter_warmup=500,                                                   # the number of warmup steps
    save_warmup=False,                                                 # we don't care about the warmup
    inits={"H0": normal(loc=70, scale=10), "Omb": normal(loc=0.3, scale=0.1), "Om": normal(loc=0.3, scale=1)},  # initial values for each parameter
    parallel_chains= 4,                                                # number of chains to run at the same time
    force_one_process_per_chain= True
)

# print very useful utilities provided by default in Stan
print("============================== fit summary ==============================")
print(fit.summary())
print("============================ fit diagnostics ============================")
print(fit.diagnose())

# show the traceplot
posterior = az.from_cmdstanpy(posterior=fit)
az.plot_trace(posterior, var_names=('H0', 'Omb', "Om"), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show corner plot
samples = [fit.stan_variable("H0"), fit.stan_variable("Omb"), fit.stan_variable("Om")]
mcsamples = MCSamples(samples=samples, names=['H0', 'Omb', "Om"], labels=['', 'Omb', "Om"])
g = plots.get_subplot_plotter()
g.triangle_plot(mcsamples, filled=True)
plt.show()
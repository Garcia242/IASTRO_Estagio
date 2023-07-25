# imports
from getdist import plots, MCSamples
from cmdstanpy import CmdStanModel
import matplotlib.pyplot as plt
from numpy.random import normal
import arviz as az

# create a model instance
model = CmdStanModel(
    stan_file="Model3/Model/cmb4.stan",                         # Stan model file location
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
    inits={"h": normal(loc=0.7, scale=1), "Omega_b": normal(loc=0.05, scale=0.1), "Omega_c": normal(loc=0.25, scale=1)},  # initial values for each parameter
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
az.plot_trace(posterior, var_names=('h', 'Omega_b', "Omega_c"), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show corner plot
samples = [fit.stan_variable("h"), fit.stan_variable("Omega_b"), fit.stan_variable("Omega_c")]
mcsamples = MCSamples(samples=samples, names=['h', 'Omega_b', "Omega_c"], labels=['h', 'Omega_b', "Omega_c"])
g = plots.get_subplot_plotter()
g.triangle_plot(mcsamples, filled=True)
plt.show()
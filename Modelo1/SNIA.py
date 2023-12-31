# imports
from getdist import plots, MCSamples
from cmdstanpy import CmdStanModel
import matplotlib.pyplot as plt
from numpy.random import normal
import arviz as az

# create a model instance
model = CmdStanModel(
    stan_file="Modelo1/Model/SNIA.stan",                         # Stan model file location
    cpp_options={
        #"STAN_NO_RANGE_CHECKS": "TRUE",  # don't check for elements out of bounds
        "STAN_THREADS": "TRUE",          # run multiple chains in parallel
        "STAN_CPP_OPTIMS": "TRUE"        # optimizations recommended by the Stan development team
    }
)
# configure and fit the model
fit = model.sample(
    data="Modelo1/Data/snia.json",                                # the location of the data file
    output_dir="output1/snia",
    iter_sampling=500,                                                 # the number of sampling steps
    iter_warmup=500,                                                   # the number of warmup steps
    save_warmup=False,                                                 # we don't care about the warmup
    inits={ "Om": normal(loc=0.3, scale=0.1), "zeta":1},  # initial values for each parameter
    parallel_chains= 4                                                # number of chains to run at the same time
)

# print very useful utilities provided by default in Stan
print("============================== fit summary ==============================")
print(fit.summary())
print("============================ fit diagnostics ============================")
print(fit.diagnose())

# show the traceplot
posterior = az.from_cmdstanpy(posterior=fit)
az.plot_trace(posterior, var_names=('Om', "zeta"), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show corner plot
samples = [ fit.stan_variable("Om"), fit.stan_variable("zeta")]
mcsamples = MCSamples(samples=samples, names=[ "Om", "zeta"], labels=[ "Om", "zeta"])
g = plots.get_subplot_plotter()
g.triangle_plot(mcsamples, filled=True)
plt.show()
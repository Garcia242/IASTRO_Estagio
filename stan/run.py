# imports
from getdist import plots, MCSamples
from cmdstanpy import CmdStanModel
import matplotlib.pyplot as plt
from numpy.random import normal
import arviz as az

# create a model instance
model = CmdStanModel(
    stan_file="Desktop/stan/model/dL.stan",                         # Stan model file location
    cpp_options={'STAN_THREADS': 'TRUE', 'STAN_CPP_OPTIMS': 'TRUE'}   # a few optimizations, feel free to ignore this
)

# configure and fit the model
fit = model.sample(
    data="Desktop/stan/data/dL.json",                                # the location of the data file
    iter_sampling=1000,                                                 # the number of sampling steps
    iter_warmup=500,                                                   # the number of warmup steps
    save_warmup=False,                                                 # we don't care about the warmup
    inits={"m": normal(loc=2, scale=1), "b": normal(loc=1, scale=1)},  # initial values for each parameter
    parallel_chains= 4                                                # number of chains to run at the same time
)

# print very useful utilities provided by default in Stan
print("============================== fit summary ==============================")
print(fit.summary())
print("============================ fit diagnostics ============================")
print(fit.diagnose())

# show the traceplot
posterior = az.from_cmdstanpy(posterior=fit)
az.plot_trace(posterior, var_names=('h', 'Om'), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show corner plot
samples = [fit.stan_variable("h"), fit.stan_variable("Om")]
mcsamples = MCSamples(samples=samples, names=["h", "Om"], labels=["h", "Om"])
g = plots.get_subplot_plotter()
g.triangle_plot(mcsamples, filled=True)
plt.show()

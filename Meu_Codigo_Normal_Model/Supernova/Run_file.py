# imports
from getdist import plots, MCSamples
from cmdstanpy import CmdStanModel
import matplotlib.pyplot as plt
from numpy.random import normal
import arviz as az

# create a model instance
model = CmdStanModel(
    stan_file="Meu_Codigo_Normal_Model/Supernova/Model/Supernova_dL.stan",                         # Stan model file location
    cpp_options={'STAN_THREADS': 'TRUE', 'STAN_CPP_OPTIMS': 'TRUE'}   # a few optimizations, feel free to ignore this
)

# configure and fit the model
fit = model.sample(
    data="Meu_Codigo_Normal_Model/Supernova/Data/snia.json",                                # the location of the data file
    output_dir="normal/snia",
    iter_sampling=500,                                                 # the number of sampling steps
    iter_warmup=500,                                                   # the number of warmup steps
    save_warmup=False,                                                 # we don't care about the warmup
    inits={"Om": normal(loc=0.3, scale=0.1)},  # initial values for each parameter
    parallel_chains= 4                                                # number of chains to run at the same time
)

# print very useful utilities provided by default in Stan
print("============================== fit summary ==============================")
print(fit.summary())
print("============================ fit diagnostics ============================")
print(fit.diagnose())

# show the traceplot
posterior = az.from_cmdstanpy(posterior=fit)
az.plot_trace(posterior, var_names=('Om'), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show corner plot
samples = [ fit.stan_variable("Om")]
mcsamples = MCSamples(samples=samples, names=["Om"], labels=["Om"])
g = plots.get_subplot_plotter()
g.triangle_plot(mcsamples, filled=True)
plt.show()
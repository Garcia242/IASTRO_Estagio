# imports
from getdist import plots, MCSamples
from cmdstanpy import CmdStanModel
import matplotlib.pyplot as plt
from numpy.random import normal
import arviz as az
from getdist import plots, gaussian_mixtures

# create a model instance
model1 = CmdStanModel(
    stan_file="Meu_Codigo/AllTogetherNow/Model/together_appart.stan",                         # Stan model file location
    cpp_options={'STAN_THREADS': 'TRUE', 'STAN_CPP_OPTIMS': 'TRUE'}   # a few optimizations, feel free to ignore this
)

# configure and fit the model
fit1 = model1.sample(
    data="Meu_Codigo/AllTogetherNow/Data/together.json",                                # the location of the data file
    iter_sampling=500,                                                 # the number of sampling steps
    iter_warmup=500,                                                   # the number of warmup steps
    save_warmup=False,                                                 # we don't care about the warmup
    inits={"H0": normal(loc=70, scale=10), "Om": normal(loc=0.3, scale=0.1)},  # initial values for each parameter
    parallel_chains= 4                                                # number of chains to run at the same time
)

model2 = CmdStanModel(
    stan_file="Meu_Codigo/Cosmic_Chronometers/Model/Chrono.stan",                         # Stan model file location
    cpp_options={'STAN_THREADS': 'TRUE', 'STAN_CPP_OPTIMS': 'TRUE'}   # a few optimizations, feel free to ignore this
)

# configure and fit the model
fit2 = model2.sample(
    data="Meu_Codigo/Cosmic_Chronometers/Data/cc.json",                                # the location of the data file
    iter_sampling=500,                                                 # the number of sampling steps
    iter_warmup=500,                                                   # the number of warmup steps
    save_warmup=False,                                                 # we don't care about the warmup
    inits={"H0": normal(loc=70, scale=10), "Om": normal(loc=0.3, scale=0.1)},  # initial values for each parameter
    parallel_chains= 4                                                # number of chains to run at the same time
)

# create a model instance
model3 = CmdStanModel(
    stan_file="Meu_Codigo/Supernova/Model/Supernova_dL.stan",                         # Stan model file location
    cpp_options={'STAN_THREADS': 'TRUE', 'STAN_CPP_OPTIMS': 'TRUE'}   # a few optimizations, feel free to ignore this
)

# configure and fit the model
fit3 = model3.sample(
    data="Meu_Codigo/Supernova/Data/snia.json",                                # the location of the data file
    iter_sampling=500,                                                 # the number of sampling steps
    iter_warmup=500,                                                   # the number of warmup steps
    save_warmup=False,                                                 # we don't care about the warmup
    inits={"Om": normal(loc=0.3, scale=0.1)},  # initial values for each parameter
    parallel_chains= 4                                                # number of chains to run at the same time
)

# print very useful utilities provided by default in Stan
print("============================== fit summary ==============================")
print(fit1.summary())
print("============================ fit diagnostics ============================")
print(fit1.diagnose())

# print very useful utilities provided by default in Stan
print("============================== fit summary ==============================")
print(fit2.summary())
print("============================ fit diagnostics ============================")
print(fit2.diagnose())

# print very useful utilities provided by default in Stan
print("============================== fit summary ==============================")
print(fit3.summary())
print("============================ fit diagnostics ============================")
print(fit3.diagnose())

# show the traceplot
posterior = az.from_cmdstanpy(posterior=fit1)
az.plot_trace(posterior, var_names=('H0', 'Om'), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show the traceplot
posterior = az.from_cmdstanpy(posterior=fit2)
az.plot_trace(posterior, var_names=('H0', 'Om'), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show the traceplot
posterior = az.from_cmdstanpy(posterior=fit3)
az.plot_trace(posterior, var_names=('Om'), compact=False, combined=False)
plt.tight_layout()
plt.show()

# show corner plot


samples1 = [fit1.stan_variable("H0"), fit1.stan_variable("Om")]
mcsamples1 = MCSamples(samples=samples1, names=["H0", "Om"], labels=["H0", "Om"])
g = plots.get_subplot_plotter()
g.triangle_plot(mcsamples1, filled=True)



# show corner plot
samples2 = [fit2.stan_variable("H0"), fit2.stan_variable("Om")]
mcsamples2 = MCSamples(samples=samples2, names=["\\H_0", "\\Omega_m"], labels=["\\H_0", "\\Omega_m"])
g = plots.get_subplot_plotter()
g.triangle_plot([mcsamples1 ,mcsamples2], filled=True)



# # show corner plot
# samples3 = [ fit3.stan_variable("Om")]
# mcsamples = MCSamples(samples=samples3, names=["Om"], labels=["Om"])
# g = plots.get_subplot_plotter()
# g.triangle_plot(mcsamples, filled=True)
# plt.show()

# samples1, samples2 = gaussian_mixtures.randomTestMCSamples(ndim=2, nMCSamples=2)
# g = plots.get_single_plotter(width_inch=4)
# g.plot_2d([samples1, samples2], ['H0','Om'], filled=False)

# samples3 = [ fit1.stan_variable("H0"), fit1.stan_variable("Om"), fit2.stan_variable("H0"), fit2.stan_variable("Om"), fit3.stan_variable("Om")]
# mcsamples = MCSamples(samples=samples3, names=["H0", "Om","H0", "Om","Om"], labels=["H0", "Om","H0", "Om", "Om"])
# g = plots.get_subplot_plotter()
# g.triangle_plot(mcsamples, filled=True)


# plt.show()
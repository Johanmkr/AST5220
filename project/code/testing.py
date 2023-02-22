import numpy as np
import os
import tk
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
import seaborn as sns
from IPython import embed
import astropy.constants as const
from astropy import units

# plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('font', family='Serif')
sns.set_theme()

# other rc parameters
plt.rc('figure', figsize=(12,7))
SMALL_SIZE = 22
MEDIUM_SIZE = 26
BIGGER_SIZE = 30
plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# folder paths
here = os.path.abspath(".")
data_path = here + "/data/"
plot_path = here + "/../output/figures/"
latex_path = here + "/../latex/"

temp_output_path = "../output/"

header_names = ["x", "eta", "Hp", "dHp", "ddHp", "OmegaB", "OmegaCDM", "OmegaLambda", "OmegaR", "OmegaNu", "OmegaK", "d_L", "dummy"]

cosmology = pd.read_csv("data/backgroundcosmology.txt", delimiter=" ", names=header_names, skiprows=1)

def testing_Omegas():
    """
    Plot for testing the Omegas evolution
    """
    xvals = cosmology["x"]
    OmegaRel = cosmology["OmegaR"] + cosmology["OmegaNu"]
    OmegaM = cosmology["OmegaB"] + cosmology["OmegaCDM"]
    omegaFig, ax1 = plt.subplots()
    ax1.plot(xvals, OmegaRel, label=r"$\Omega_\mathrm{rel}$", color="red")
    ax1.plot(xvals, OmegaM, label=r"$\Omega_\mathrm{M}$", color="green")
    ax1.plot(xvals, cosmology["OmegaLambda"], label=r"$\Omega_\Lambda$", color="orange")
    ax1.plot(xvals, cosmology["OmegaLambda"]+OmegaRel+OmegaM, label="Sum", color="blue")
    ax1.set_xlabel(r"$x$")
    ax1.set_title(r"Density fractions $\Omega_X$", loc="left")
    omegaFig.legend(loc='center right', bbox_to_anchor=(.95, 0.5),
          ncol=1, fancybox=True)
    omegaFig.tight_layout()
    omegaFig.savefig(plot_path+"omegas.pdf", bbox_inches=None)
    # omegaFig.show()

def testing_Hp():
    """
    Plot for testing Hp
    """
    xvals = cosmology["x"]
    Hp = cosmology["Hp"]
    dHp = cosmology["dHp"]
    ratio = dHp/Hp 
    HpFig, ax1 = plt.subplots()
    ax1.plot(xvals, ratio, color="blue", label=r"$\frac{1}{\mathcal{H}}\frac{\mathrm{d}\mathcal{H}}{\mathrm{d}x}$")
    ax1.set_title(r"First sanity check of $\mathcal{H}(x)$", loc="left")
    ax1.set_xlabel(r"$x$")
    ax1.legend(loc="upper left", fancybox=True)
    HpFig.tight_layout()
    HpFig.savefig(plot_path + "Hp_test.pdf", bbox_inches=None)

def testing_Hp2():
    """
    Plot for testing Hp (second)
    """
    xvals = cosmology["x"]
    Hp = cosmology["Hp"]
    ddHp = cosmology["ddHp"]
    ratio = ddHp/Hp 
    HpFig, ax1 = plt.subplots()
    ax1.plot(xvals, ratio, color="blue", label=r"$\frac{1}{\mathcal{H}}\frac{\mathrm{d}^2\mathcal{H}}{\mathrm{d}x^2}$")
    ax1.set_title(r"Second sanity check of $\mathcal{H}(x)$", loc="left")
    ax1.set_xlabel(r"$x$")
    ax1.legend(loc="lower left", fancybox=True)
    HpFig.tight_layout()
    HpFig.savefig(plot_path + "Hp_test2.pdf", bbox_inches=None)



if __name__=="__main__":
    # testing_Omegas()
    testing_Hp()
    testing_Hp2()
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
sns.color_palette("hls", 8)

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

header_names = ["x", "eta", "deta", "Hp", "dHp", "ddHp", "OmegaB", "OmegaCDM", "OmegaLambda", "OmegaR", "OmegaNu", "OmegaK", "d_L", "dummy"]

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
    ddHp = cosmology["ddHp"]
    ratiodHp = dHp/Hp 
    ratioddHp = ddHp/Hp 
    HpFig, ax1 = plt.subplots()
    ax1.plot(xvals, ratioddHp, color="red", label=r"$\frac{1}{\mathcal{H}}\frac{\mathrm{d}^2\mathcal{H}}{\mathrm{d}x^2}$")
    ax1.plot(xvals, ratiodHp, color="blue", label=r"$\frac{1}{\mathcal{H}}\frac{\mathrm{d}\mathcal{H}}{\mathrm{d}x}$")
    ax1.set_title(r"Sanity check of $\mathcal{H}(x)$", loc="left")
    ax1.set_xlabel(r"$x$")
    ax1.legend(loc="center left", fancybox=True)
    HpFig.tight_layout()
    HpFig.savefig(plot_path + "Hp_test.pdf", bbox_inches=None)

def testing_eta():
    xvals = cosmology["x"].where(cosmology["x"]<0)
    eta = cosmology["eta"].where(cosmology["x"]<0)
    deta = cosmology["deta"].where(cosmology["x"]<0)
    Hp = cosmology["Hp"].where(cosmology["x"]<0)
    # xvals = cosmology["x"]
    # eta = cosmology["eta"]
    # deta = cosmology["deta"]
    # Hp = cosmology["Hp"]
    c = const.c

    etaHpc = eta*Hp/c 
    detaHpc = deta*Hp/c

    etaFig, ax1 = plt.subplots()
    ax1.plot(xvals, etaHpc, color="blue", label=r"$\frac{\eta\mathcal{H}}{c}$")
    ax1.plot(xvals, detaHpc, color="red", label=r"$\frac{\mathcal{H}}{c}\frac{\mathrm{d}\eta}{\mathrm{d}x}$")
    ax1.set_title(r"Sanity check for $\eta(x)$", loc="left")
    ax1.set_xlabel(r"$x$")
    # ax1.set_yscale("log")
    ax1.legend(loc="best", fancybox=True)
    etaFig.tight_layout()
    etaFig.savefig(plot_path + "eta_test.pdf", bbox_inches=None)



if __name__=="__main__":
    # testing_Omegas()
    # testing_Hp()
    testing_eta()
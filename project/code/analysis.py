import numpy as np
import os
import tkinter
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
import seaborn as sns
from IPython import embed
import astropy.constants as const
from astropy import units
os.environ["QT_QPA_PLATFORM"] = "wayland"

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

header_names = ["x", "eta", "Hp", "dHp", "OmegaB", "OmegaCDM", "OmegaLambda", "OmegaR", "OmegaNu", "OmegaK", "d_L", "dummy"]

cosmology = pd.read_csv("data/backgroundcosmology.txt", delimiter=" ", names=header_names)


# sfit = pd.read_csv(data_path + "results_supernovafitting.txt", delimiter=" ", names=["chi2", "h", "OmegaM", "OmegaK"])

##############
#   Omega plot
##############

# Omega_gamma = 2*np.pi**2 / 30 * (const.k_b*2.7255)**4 / const.hbar**3 / const.c*5 * 8*np.pi * const.G /

def omega_restrictions_plot():
    sfit = np.loadtxt(data_path+"results_supernovafitting.txt", skiprows=300)
    chi2, h, OmegaM, OmegaK = sfit[:,0], sfit[:,1], sfit[:,2], sfit[:,3]

    oneSDthres = 3.53
    chi2min = np.min(chi2)

    accepted = (chi2 - chi2min) < oneSDthres

    selected_omegaM = np.where(accepted, OmegaM, np.nan)
    selected_omegaK = np.where(accepted, OmegaK, np.nan)
    selected_omegaM = selected_omegaM[np.isfinite(selected_omegaM)]
    selected_omegaK = selected_omegaK[np.isfinite(selected_omegaK)]


    fig, ax = plt.subplots()

    ax.scatter(selected_omegaM, 1-(selected_omegaK+selected_omegaM))
    fig.show()
    # fig.savefig(plot_path+"omega_restrictions.pdf", bbox_inches=None)



# print(sfit)



if __name__=="__main__":
    omega_restrictions_plot()
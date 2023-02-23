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

# Visual parameters for saving/showing
SAVE = True
PUSH = False
SHOW = False 
TIGHT = True

# folder paths
here = os.path.abspath(".")
data_path = here + "/data/"
plot_path = here + "/../output/figures/"
latex_path = here + "/../latex/"

temp_output_path = "../output/"





#   Dictionary with labels
lbls = {
    "x": r"$x$",
    "eta": r"$\eta$",
    "deta": r"$\frac{\mathrm{d}\eta}{\mathrm{d}x}$",
    "Hp": r"$\mathcal{H}$",
    "dHp": r"$\frac{\mathrm{d}\mathcal{H}}{\mathcal{d}x}$",
    "ddHp": r"$\frac{\mathrm{d}^2\mathcal{H}}{\mathcal{d}x^2}$",
    "OmegaB": r"$\Omega_b$",
    "OmegaCDM": r"$\Omega_\mathrm{CDM}$",
    "OmegaLambda": r"$\Omega_\Lambda$",
    "OmegaR": r"$\Omega_\gamma$",
    "OmegaNu": r"$\Omega_\nu$",
    "OmegaK": r"$\Omega_k$",
    "OmegaRad": r"$\Omega_mathrm{rad}$",
    "OmegaM": r"$\Omega_\mathrm{M}$",
    "d_L": r"$d_L$",
    "chi2": r"$\chi^2$",
    "H": r"$H$",
    "h": r"$h$",
    "t": r"$t$",
    "z": r"$z$"
}

#   Set up data class
class Data:
    def __init__(self, filename, skiprows=0):
        f = open(data_path+filename)
        header = f.readline()
        f.close()
        column_names = header.split()[1:]
        # Accound for the units in supernovadata.txt
        for i, name in reversed(list(enumerate(column_names))):
            if name in ["(Gpc)", "Acceptrate"]:
                column_names.pop(i)
        data = np.loadtxt(data_path + filename, skiprows=skiprows)
        data_dict = {}
        for i, key in enumerate(column_names):
            data_dict[key] = data[:,i]
        self.dF = pd.DataFrame(data_dict)
    
    def __call__(self):
        return self.dF
    
    def __str__(self):
        return self.to_string()
    
    def __getitem__(self, item):
        return self.dF[item]
    
    def print_frame(self):
        print(self.dF)


#   Generate data objects

Cosmology = Data("backgroundcosmology.txt")
Sdata = Data("supernovadata.txt")
Sfit = Data("results_supernovafitting.txt", skiprows=300)

def save_push(fig, pdf_name, save=SAVE, push=PUSH, show=SHOW, tight=TIGHT):
    if tight:
        fig.tight_layout()
    file = plot_path + pdf_name.replace('.pdf', '').strip() + ".pdf"
    if save and pdf_name != 'none':
        print(f'Saving plot: {file}')
        fig.savefig(file, bbox_inches="tight")
    if push and pdf_name != 'none':
        os.system(f"git add {file}")
        os.system("git commit -m 'upload plot'")
        os.system("git push")
    if show:
        plt.show()
    else:
        plt.clf()

    plt.close()


def conformal_hubble_factor():
    xvals = Cosmology["x"]
    Hp = Cosmology["Hp"]
    chf, ax = plt.subplots()
    ax.plot(xvals, Hp, color="blue", label=lbls["Hp"])
    ax.set_xlabel(lbls["x"])
    ax.set_ylabel("insert unit")
    ax.set_title(r"Conformal Hubble factor $\mathcal{H}(x)$", loc="left")
    ax.legend(loc="best", fancybox=True)
    ax.set_yscale("log")
    save_push(chf, "conformal_hubble_factor")


def cosmic_time():
    xvals = Cosmology["x"]
    t = Cosmology["t"] * units.s.to("Gyr")
    ct, ax = plt.subplots()
    ax.plot(xvals, t, color="blue", label=lbls["t"])
    ax.set_xlabel(lbls["x"])
    ax.set_ylabel(lbls["t"]+" [Gyr]")
    ax.set_title(r"Cosmic time $t(x)$", loc="left")
    ax.legend(loc="best", fancybox=True)
    ax.set_yscale("log")
    save_push(ct, "cosmic_time")

def supernova_data():
    zvals_sn = Sdata["z"]
    dL_obs = Sdata["d_L"]
    error = Sdata["Error"]
    xvals_pred = Cosmology["x"]
    dL_pred = Cosmology["d_L"]*units.m.to("Gpc")
    zvals_pred = np.exp(-xvals_pred)-1

    sdFig, ax = plt.subplots()
    ax.errorbar(zvals_sn, dL_obs, yerr=error, label="Observation", fmt="none", ecolor="red", marker="s", ms=20, mew=5, capsize=.5, lw=1)
    ax.plot(zvals_pred, dL_pred, color="blue", label="Prediction")
    ax.set_yscale("log")
    ax.set_ylim(0.1,12)
    ax.set_xlim(0,1.4)
    ax.set_xlabel(lbls["z"])
    ax.set_ylabel(lbls["d_L"] + " [Gpc]")
    ax.set_title(r"Luminosity distance $d_L$", loc="left")
    ax.legend(loc="best", fancybox=True)    
    
    save_push(sdFig, "supernova_data")



def omega_restrictions_plot():
    chi2, h, OmegaM, OmegaK = sfit[:,0], sfit[:,1], sfit[:,2], sfit[:,3]
    chi2 

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



if __name__=="__main__":
    # conformal_hubble_factor()
    # cosmic_time()
    supernova_data()
    # omega_restrictions_plot()
    # embed()

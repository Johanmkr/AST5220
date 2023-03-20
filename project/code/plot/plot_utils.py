import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython import embed
import astropy.constants as const
from astropy import units


#   Colouring for cosmology project
Colors = {
    "OmegaRad": "orange",
    "OmegaM": "green",
    "OmegaLambda": "purple",
    "etaHp/c":   "blue",
    "Hp/cdetadx":    "red",
    "ddHpddx":  "blue",
    "dHpdx": "red",
    "Hp": "blue",
    "t": "blue",
    "eta/c": "red",
    "d_L_obs": "red",
    "d_L_fid": "blue",
    "d_L_best": "springgreen",
    "hist": "turquoise",
    "gaussian": "forestgreen"
}
# colors = plt.cm.winter(np.linspace(0,1,5))
CMAP = "winter"
 
# plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
sns.set_theme(palette="winter")

# other rc parameters
plt.rc('figure', figsize=(12,7))
SMALL_SIZE = 30
MEDIUM_SIZE = 35
BIGGER_SIZE = 40
plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('axes', grid=False)
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('lines', linewidth=3)

# Visual parameters for saving/showing
SAVE = False
PUSH = False
SHOW = True
TIGHT = True

# folder paths
here = os.path.abspath(".")
data_path = here + "/../data/"
plot_path = here + "/../../output/figures/"
latex_path = here + "/../../tex/tables/"

def save_push(fig, pdf_name, save=SAVE, push=PUSH, show=SHOW, tight=TIGHT):
    """Function to save, plot and/push figures

    Args:
        fig (plt.figure): Figure object.
        pdf_name (str): Name of pdf
        save (bool, optional): Save? Defaults to SAVE.
        push (bool, optional): Push to github? Defaults to PUSH.
        show (bool, optional): Show plot? Defaults to SHOW.
        tight (bool, optional): Tight layout on plot? Defaults to TIGHT.
    """
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
        if filename[-3:] == "txt":
            f = open(data_path+filename)
            header = f.readline()[1:]
            f.close()
            column_names = header.split(" ")
            actual_column_names = []
            for name in column_names:
                if name == "":
                    pass
                elif name[0] in ["(", "["]:
                    pass
                else:
                    actual_column_names.append(name)
            data = np.loadtxt(data_path + filename, skiprows=skiprows)
            data_dict = {}
            for i, key in enumerate(actual_column_names):
                data_dict[key] = data[:,i]
            self.dF = pd.DataFrame(data_dict)
        elif filename[-3:] == "csv":
            self.dF = pd.read_csv(data_path+filename, skiprows=range(1,skiprows+1))
            self.dF.columns = self.dF.columns.str.strip()

        else:
            print("----Provide valid file!----")
    
    def __call__(self):
        return self.dF
    
    def __str__(self):
        return self.dF.to_string()
    
    def __getitem__(self, item):
        try:
            return self.dF[item]
        except KeyError:
            key_vals = self.dF.keys()
            for key in key_vals:
                if item in key:
                    return self.dF[key]
    
    def print_frame(self):
        print(self.dF)

#   Different utility functions

# def add_shaded_areas(axis):
    
if __name__=="__main__":
    print("Utility file only.")
    
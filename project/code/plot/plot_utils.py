import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
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
SMALL_SIZE = 25
MEDIUM_SIZE = 30
BIGGER_SIZE = 35
plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('lines', linewidth=2)

# Visual parameters for saving/showing
SAVE = True
PUSH = False
SHOW = False 
TIGHT = True

# misc
CMAP = "winter"

# folder paths
here = os.path.abspath(".")
data_path = here + "/../data/"
plot_path = here + "/../../output/figures/"
latex_path = here + "/../../latex/"

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

if __name__=="__main__":
    print("Utility file only.")
    
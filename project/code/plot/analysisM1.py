from plot_utils import *

#   Generate data objects
Cosmology = Data("backgroundcosmology.csv")
Sdata = Data("supernovadata.txt")
Sfit = Data("results_supernovafitting.csv", skiprows=300)
ValueTab = Data("table_of_values.csv")
# embed()




def conformal_hubble_factor():
    """Plot the conformal Hubble factor Hp agains x.
    """
    xvals = Cosmology["x"]
    Hp = Cosmology["Hp"]

    chf, ax = plt.subplots()
    ax.plot(xvals, Hp, color=colors[0], label=lbls["Hp"])
    ax.set_xlabel(lbls["x"])
    ax.set_ylabel(r"$\mathcal{H}$ [s$^{-1}$]")
    ax.set_title(r"Conformal Hubble factor $\mathcal{H}(x)$", loc="left")
    # ax.legend(loc="best", fancybox=True)
    ax.set_yscale("log")

    save_push(chf, "conformal_hubble_factor")

def cosmic_conformal_time():
    """Make a plot of the cosmic time t against x.
    """
    xvals = Cosmology["x"]
    t = Cosmology["t"] * units.s.to("Gyr")

    eta = Cosmology["eta"]*units.m
    eta_c = eta/const.c.to("m/Gyr")

    ct, ax = plt.subplots()
    ax.plot(xvals, t, color=colors[0], label=lbls["t"])
    ax.plot(xvals, eta_c, color=colors[-1], label=r"$\frac{\eta}{c}$")

    ax.set_xlabel(lbls["x"])
    ax.set_ylabel(lbls["t"]+" [Gyr]")
    ax.set_title(r"Cosmic time $t(x)$ and conformal time $\eta(x)/c$.", loc="left")
    ax.legend(loc="best", fancybox=True)
    # ax.set_yscale("log")
    save_push(ct, "cosmic_conformal_time")


def supernova_data():
    """Plot predicted luminosity distance together with the observed supernova data.
    """
    zvals_sn = Sdata["z"]
    dL_obs = Sdata["d_L"]
    error = Sdata["Error"]
    xvals_pred = Cosmology["x"]
    dL_pred = Cosmology["d_L"]*units.m.to("Gpc")
    zvals_pred = np.exp(-xvals_pred)-1

    sdFig, ax = plt.subplots()
    ax.errorbar(zvals_sn, dL_obs/zvals_sn, yerr=error/zvals_sn, label="Observation", fmt="none",  ecolor=colors[0])
    ax.plot(zvals_pred, dL_pred/zvals_pred, color=colors[-1], label="Prediction")
    # ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylim(3.5,8)
    ax.set_xlim(0.005,1.45)
    ax.set_xlabel(lbls["z"])
    ax.set_ylabel(r"$d_L/z$" + " [Gpc]")
    ax.set_title(r"Luminosity distance $d_L$", loc="left")
    ax.legend(loc="upper left", fancybox=True)    
    
    save_push(sdFig, "supernova_data")



def omega_restrictions_plot():
    """Make a plot of the 1SD confidence of the chi2 values in the OmegaM, OmegaLambda plane.
    """
    chi2 = Sfit["chi2"]
    OmegaM = Sfit["OmegaM"]
    OmegaK = Sfit["OmegaK"]

    oneSDthres = 3.53
    chi2min = np.min(chi2)

    accepted = (chi2 - chi2min) < oneSDthres

    selected_omegaM = np.where(accepted, OmegaM, np.nan)
    selected_omegaK = np.where(accepted, OmegaK, np.nan)
    selected_chi2 = np.where(accepted, chi2, np.nan)
    selected_omegaM = selected_omegaM[np.isfinite(selected_omegaM)]
    selected_omegaK = selected_omegaK[np.isfinite(selected_omegaK)]
    selected_chi2 = selected_chi2[np.isfinite(selected_chi2)]
    selected_omegaLambda = 1 - (selected_omegaK + selected_omegaM)

    OmegaPlane, ax = plt.subplots() 
    scat = ax.scatter(selected_omegaM, selected_omegaLambda, c=selected_chi2, cmap=CMAP)
    flatline = ax.plot((0,1), (1,0), color="black", ls="--")
    # l1 = ax.legend(loc="upper right", fancybox=True, bbox_to_anchor=(1, 1.15))

    cbar = plt.colorbar(scat, ax=ax, label=r"$\chi^2$")

    ax.set_xlabel(lbls["OmegaM"])
    ax.set_ylabel(lbls["OmegaLambda"])
    ax.set_ylim(.2,1)
    ax.set_xlim(0,0.6)
    ax.set_title(r"$1\sigma$ confidence plot of $(\Omega_M \;\mathrm{x}\; \Omega_\Lambda)$", loc="left")
    # fig.savefig(plot_path+"omega_restrictions.pdf", bbox_inches=None)

    save_push(OmegaPlane, "omega_plane")

def posterior_pdf():
    h = Sfit["h"]

    #   Make H0 from h
    H0 = 100* h

    #create count and bins
    counts, bins = np.histogram(H0, bins=150)

    # Gaussian func
    sigma, mu = np.std(H0), np.mean(H0)
    gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))

    ppdf, ax = plt.subplots()
    hstg = ax.hist(bins[:-1], bins, weights=counts, color=colors[-1], label="Samples", density=True)
    ax.plot(bins, gaussian, color=colors[0], label="Fitted pdf")
    ax.set_xlabel(r"$H_0$ [km s$^{-1}$Mpc$^{-1}$]")
    ax.set_ylabel("Probability")
    ax.set_title(r"Posterior pdf of $H_0$", loc="left")
    l1 = ax.legend(loc="best", fancybox=True)

    save_push(ppdf, "posterior_pdf")

def create_table():
    styler = ValueTab.dF.style
    styler.format({"x": '{:.2f}', "z": '{:.2f}', "t [Gyr]": '{:.6f}'})
    styler.hide(axis="index")
    styler.to_latex(latex_path + "value_table.tex")



if __name__=="__main__":
    conformal_hubble_factor()
    cosmic_conformal_time()
    supernova_data()
    omega_restrictions_plot()
    posterior_pdf()
    create_table()

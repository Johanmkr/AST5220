from plot_utils import *

#   Generate data objects
Cosmology = Data("backgroundcosmology.csv")
Sdata = Data("supernovadata.txt")
Sfit = Data("results_supernovafitting.csv", skiprows=300)
ValueTab = Data("table_of_values.csv")

x_min = Cosmology["x"][0]
x_max = np.asarray(Cosmology["x"])[-1]
x_RM = ValueTab["x"][0]
x_ML = ValueTab["x"][1]
tol = .5e-1

def set_regimes(ax):
    rad_area = ax.axvspan(x_min, x_RM-tol, color=Colors["OmegaRad"], alpha=.1, label=r"$\Omega_\mathrm{rad}$")
    mat_area = ax.axvspan(x_RM+tol, x_ML-tol, color=Colors["OmegaM"], alpha=.1, label=r"$\Omega_\mathrm{M}$")
    lam_area = ax.axvspan(x_ML+tol, x_max, color=Colors["OmegaLambda"], alpha=.1, label=r"$\Omega_\Lambda$")
    return [rad_area, mat_area, lam_area]

def testing_Omegas():
    """
    Plot for testing the Omegas evolution
    """
    xvals = Cosmology["x"]
    OmegaRad = Cosmology["OmegaR"] + Cosmology["OmegaNu"]
    OmegaM = Cosmology["OmegaB"] + Cosmology["OmegaCDM"]
    
    omegaFig, ax1 = plt.subplots()
    ax1.plot(xvals, Cosmology["OmegaLambda"]+OmegaRad+OmegaM, label="Sum", ls="--", color="black")
    ax1.plot(xvals, OmegaRad, label=r"$\Omega_\mathrm{rad}$", color=Colors["OmegaRad"])
    ax1.plot(xvals, OmegaM, label=r"$\Omega_\mathrm{M}$", color=Colors["OmegaM"])
    ax1.plot(xvals, Cosmology["OmegaLambda"], label=r"$\Omega_\Lambda$", color=Colors["OmegaLambda"])
    ax1.set_xlabel(r"$x$")
    ax1.set_title(r"Density fractions $\Omega_X$", loc="left")
    ax1.legend(loc='center left', ncol=1, fancybox=True)
    
    save_push(omegaFig, "testing_omegas")

def testing_Hp():
    """
    Plot for testing Hp
    """

    xvals = Cosmology["x"]
    Hp = Cosmology["Hp"]
    dHp = Cosmology["dHp"]
    ddHp = Cosmology["ddHp"]
    ratiodHp = dHp/Hp 
    ratioddHp = ddHp/Hp 

    HpFig, ax1 = plt.subplots()
    line1, = ax1.plot(xvals, ratioddHp, label=r"$\frac{1}{\mathcal{H}}\frac{\mathrm{d}^2\mathcal{H}}{\mathrm{d}x^2}$", color=Colors["ddHpddx"])
    line2, = ax1.plot(xvals, ratiodHp, label=r"$\frac{1}{\mathcal{H}}\frac{\mathrm{d}\mathcal{H}}{\mathrm{d}x}$", color=Colors["dHpdx"])

    regimes = set_regimes(ax1)
    
    ax1.set_title(r"Sanity check of $\mathcal{H}(x)$", loc="left")
    ax1.set_xlabel(r"$x$")
    legend1 = ax1.legend([line1, line2], [line1.get_label(), line2.get_label()], loc="center left", fancybox=True)

    # legend2 = HpFig.legend([rad_area, mat_area, lam_area], [rad_area.get_label(), mat_area.get_label(), lam_area.get_label()], loc="upper right", fancybox=True, ncol=3, bbox_to_anchor=[0.97,0.965], fontsize=24)
    legend2 = HpFig.legend(regimes, [regime.get_label() for regime in regimes], loc="upper right", fancybox=True, ncol=3, bbox_to_anchor=[0.97,0.965], fontsize=24)


    save_push(HpFig, "Hp_test")

def testing_eta():
    xvals = Cosmology["x"].where(Cosmology["x"]<0)
    eta = Cosmology["eta"].where(Cosmology["x"]<0)
    deta = Cosmology["deta"].where(Cosmology["x"]<0)
    Hp = Cosmology["Hp"].where(Cosmology["x"]<0)
    c = const.c
    etaHpc = eta*Hp/c 
    detaHpc = deta*Hp/c

    etaFig, ax1 = plt.subplots()
    line1, = ax1.plot(xvals, etaHpc, label=r"$\frac{\eta\mathcal{H}}{c}$", color=Colors["etaHp/c"])
    line2, = ax1.plot(xvals, detaHpc, label=r"$\frac{\mathcal{H}}{c}\frac{\mathrm{d}\eta}{\mathrm{d}x}$", color=Colors["Hp/cdetadx"])

    ax1.set_title(r"Sanity check for $\eta(x)$", loc="left")
    ax1.set_xlabel(r"$x$")

    rad_area = ax1.axvspan(x_min, x_RM-tol, color=Colors["OmegaRad"], alpha=.1, label=r"$\Omega_\mathrm{rad}$")
    mat_area = ax1.axvspan(x_RM+tol, x_ML-tol, color=Colors["OmegaM"], alpha=.1, label=r"$\Omega_\mathrm{M}$")
    lam_area = ax1.axvspan(x_ML+tol, x_max, color=Colors["OmegaLambda"], alpha=.1, label=r"$\Omega_\Lambda$")

    # ax1.set_yscale("log")
    # ax1.legend(loc="best", fancybox=True)
    legend1 = ax1.legend([line1, line2], [line1.get_label(), line2.get_label()], loc="upper left", fancybox=True)
    legend2 = etaFig.legend([rad_area, mat_area, lam_area], [rad_area.get_label(), mat_area.get_label(), lam_area.get_label()], loc="upper right", fancybox=True, ncol=3, bbox_to_anchor=[0.97,0.965], fontsize=24)


    save_push(etaFig, "eta_test")



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

    # Set regimes
    rad_area = ax1.axvspan(x_min, x_RM-tol, color=Colors["OmegaRad"], alpha=.1, label=r"$\Omega_\mathrm{rad}$")
    mat_area = ax1.axvspan(x_RM+tol, x_ML-tol, color=Colors["OmegaM"], alpha=.1, label=r"$\Omega_\mathrm{M}$")
    lam_area = ax1.axvspan(x_ML+tol, x_max, color=Colors["OmegaLambda"], alpha=.1, label=r"$\Omega_\Lambda$")

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
    testing_Omegas()
    testing_Hp()
    testing_eta()
    conformal_hubble_factor()
    cosmic_conformal_time()
    supernova_data()
    omega_restrictions_plot()
    posterior_pdf()
    # create_table()

from plot_utils import *

#   Generate data objects
Cosmology = Data("backgroundcosmology.csv")
Sdata = Data("supernovadata.txt")
Sfit = Data("results_supernovafitting.csv", skiprows=300)
ValueTab = Data("table_of_values.csv")

LumDistData = Data("backgroundcosmologyLumDist.csv")
LumDistDataBestFit = Data("bestFitBackground.csv")

idx_today = np.argmin(np.abs(Cosmology["x"]))

x_min = Cosmology["x"][0]
x_max = np.asarray(Cosmology["x"])[-1]
x_RM = ValueTab["x"][0]
x_ML = ValueTab["x"][1]
x_accel_start = ValueTab["x"][2]
tol = 1e-2

OmegaRad0 = Cosmology["OmegaR"][idx_today] + Cosmology["OmegaNu"][idx_today]
OmegaM0 = Cosmology["OmegaB"][idx_today] + Cosmology["OmegaCDM"][idx_today]
OmegaLambda0 = Cosmology["OmegaLambda"][idx_today]
H0 = Cosmology["Hp"][idx_today]

# embed()


def set_regimes(ax, borders=True):
    rad_area = ax.axvspan(x_min, x_RM-tol, color=Colors["OmegaRad"], alpha=.1, label=r"$\Omega_\mathrm{rad}$")
    mat_area = ax.axvspan(x_RM+tol, x_ML-tol, color=Colors["OmegaM"], alpha=.1, label=r"$\Omega_\mathrm{M}$")
    lam_area = ax.axvspan(x_ML+tol, x_max, color=Colors["OmegaLambda"], alpha=.1, label=r"$\Omega_\Lambda$")
    if borders:
        ax.axvline(x_RM, color="black", ls="--")
        ax.axvline(x_ML, color="black", ls="--")
    return [rad_area, mat_area, lam_area]

def testing_Omegas():
    """
    Plot for testing the Omegas evolution
    """
    xvals = Cosmology["x"]
    OmegaRad = Cosmology["OmegaR"] + Cosmology["OmegaNu"]
    OmegaM = Cosmology["OmegaB"] + Cosmology["OmegaCDM"]
    
    omegaFig, ax1 = plt.subplots()
    ax1.plot(xvals, OmegaRad, label=r"$\Omega_\mathrm{rad}$", color=Colors["OmegaRad"])
    ax1.plot(xvals, OmegaM, label=r"$\Omega_\mathrm{M}$", color=Colors["OmegaM"])
    ax1.plot(xvals, Cosmology["OmegaLambda"], label=r"$\Omega_\Lambda$", color=Colors["OmegaLambda"])
    ax1.plot(xvals, Cosmology["OmegaLambda"]+OmegaRad+OmegaM, label=r"$\mathrm{Sum}$", ls="--", color=Colors["analytical"], lw=2)
    ax1.set_xlabel(r"$x$")
    ax1.set_title(r"$\mathrm{Density\ fractions,}\ \Omega_i$", loc="left")
    ax1.legend(loc='center left', ncol=1, fancybox=True)
    ax1.minorticks_on()

    
    regimes = set_regimes(ax1)
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

    regimes = set_regimes(ax1, borders=False)

    #   adding the analytical values
    #   for double deriv
    ax1.hlines(1, x_min, x_RM, color=Colors["analytical"], ls="--", lw=2)
    ax1.hlines(1/4, x_RM, x_ML, color=Colors["analytical"], ls="--", lw=2)
    ax1.hlines(1, x_ML, x_max, color=Colors["analytical"], ls="--", lw=2)

    #   for single deriv
    ax1.hlines(-1, x_min, x_RM, color=Colors["analytical"], ls="--", lw=2)
    ax1.hlines(-1/2, x_RM, x_ML, color=Colors["analytical"], ls="--", lw=2)
    # ax1.hlines(1, x_ML, x_max, color="white", ls="--")
    
    ax1.set_title(r"$\mathrm{Sanity\ check\ of}\ \mathcal{H}(x)$", loc="left")
    ax1.set_xlabel(r"$x$")
    legend1 = ax1.legend([line1, line2], [line1.get_label(), line2.get_label()], loc="center left", fancybox=True)
    ax1.minorticks_on()


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

    ax1.set_title(r"$\mathrm{Sanity\ check\ for}\ \eta(x)$", loc="left")
    ax1.set_xlabel(r"$x$")

    regimes = set_regimes(ax1, borders=False)

    #   Adding analytical solutions
    ax1.hlines(1, x_min, x_RM, color=Colors["analytical"], ls="--", lw=2)
    ax1.hlines(2, x_RM, x_ML, color=Colors["analytical"], ls="--", lw=2)
    # ax1.axvline(x_ML, color="snow", lw=2, ls="--")

    rad_area = ax1.axvspan(x_min, x_RM-tol, color=Colors["OmegaRad"], alpha=.1, label=r"$\Omega_\mathrm{rad}$")
    mat_area = ax1.axvspan(x_RM+tol, x_ML-tol, color=Colors["OmegaM"], alpha=.1, label=r"$\Omega_\mathrm{M}$")
    lam_area = ax1.axvspan(x_ML+tol, x_max, color=Colors["OmegaLambda"], alpha=.1, label=r"$\Omega_\Lambda$")

    # ax1.set_yscale("log")
    # ax1.legend(loc="best", fancybox=True)
    legend1 = ax1.legend([line1, line2], [line1.get_label(), line2.get_label()], loc="upper left", fancybox=True)
    legend2 = etaFig.legend(regimes, [regime.get_label() for regime in regimes], loc="upper right", fancybox=True, ncol=3, bbox_to_anchor=[0.97,0.965], fontsize=24)


    ax1.minorticks_on()

    save_push(etaFig, "eta_test")



def conformal_hubble_factor():
    """Plot the conformal Hubble factor Hp agains x.
    """
    xvals = Cosmology["x"]
    Hp = Cosmology["Hp"]*100*(1/units.s).to("km/s/Mpc")
    xvals_rad = np.linspace(x_min, x_RM, 100)
    xvals_mat = np.linspace(x_RM, x_ML, 100)
    xvals_lam = np.linspace(x_ML, x_max, 100)
    rad_anal = H0*np.sqrt(OmegaRad0)*np.exp(-xvals_rad)*100*(1/units.s).to("km/s/Mpc")
    mat_anal = H0*np.sqrt(OmegaM0)*np.exp(-xvals_mat/2)*100*(1/units.s).to("km/s/Mpc")
    lam_anal = H0*np.sqrt(OmegaLambda0)*np.exp(xvals_lam)*100*(1/units.s).to("km/s/Mpc")

    chf, ax = plt.subplots()
    ax.plot(xvals, Hp, color=Colors["Hp"], label=lbls["Hp"])
    ax.axvline(x_accel_start, color="black", ls="--", label=r"$\mathrm{Accel.\ start}$")

    #   Plot analytical solutions in regimes
    ax.plot(xvals_rad, rad_anal, color=Colors["analytical"], ls="--", lw=2)
    ax.plot(xvals_mat, mat_anal, color=Colors["analytical"], ls="--", lw=2)
    ax.plot(xvals_lam, lam_anal, color=Colors["analytical"], ls="--", lw=2)


    ax.set_xlabel(lbls["x"])
    ax.set_ylabel(r"$100\ \mathrm{kms}^{-1}\mathrm{Mpc}^{-1}$")
    ax.set_title(r"$\mathrm{Conformal\ Hubble\ factor}\ \mathcal{H}(x)$", loc="left")
    ax.legend(loc="best", fancybox=True)
    ax.set_yscale("log")

    regimes = set_regimes(ax, borders=False)

    # Set regimes
    # rad_area = ax1.axvspan(x_min, x_RM-tol, color=Colors["OmegaRad"], alpha=.1, label=r"$\Omega_\mathrm{rad}$")
    # mat_area = ax1.axvspan(x_RM+tol, x_ML-tol, color=Colors["OmegaM"], alpha=.1, label=r"$\Omega_\mathrm{M}$")
    # lam_area = ax1.axvspan(x_ML+tol, x_max, color=Colors["OmegaLambda"], alpha=.1, label=r"$\Omega_\Lambda$")
    ax.minorticks_on()

    save_push(chf, "conformal_hubble_factor")

def cosmic_conformal_time():
    """Make a plot of the cosmic time t against x.
    """
    xvals = Cosmology["x"]
    t = Cosmology["t"] * units.s.to("Gyr")

    eta = Cosmology["eta"]*units.m
    eta_c = eta/const.c.to("m/Gyr")

    ct, ax = plt.subplots()
    ax.plot(xvals, t, color=Colors["t"], label=lbls["t"])
    ax.plot(xvals, eta_c, color=Colors["eta/c"], label=r"$\frac{\eta}{c}$")


    ax.set_xlabel(lbls["x"])
    ax.set_ylabel(r"$\mathrm{Gyr}$")
    ax.set_title(r"$\mathrm{Cosmic\ time}\ t(x)\  \mathrm{and\ conformal\ time}\ \eta(x)/c$.", loc="left")
    ax.legend(loc="best", fancybox=True)
    regimes = set_regimes(ax, borders=False)
    ax.set_yscale("log")
    ax.minorticks_on()

    save_push(ct, "cosmic_conformal_time")


def supernova_data():
    """Plot predicted luminosity distance together with the observed supernova data.
    """
    zvals_sn = Sdata["z"]
    dL_obs = Sdata["d_L"]
    error = Sdata["Error"]

    #   Fiducial cosmology
    xvals_pred = LumDistData["x"]
    dL_pred = LumDistData["d_L"]*units.m.to("Gpc")
    zvals_pred = np.exp(-xvals_pred)-1

    #   Best fit cosmology
    xvalsBF = LumDistDataBestFit["x"]
    dLBF = LumDistDataBestFit["d_L"]*units.m.to("Gpc")
    zvalsBF = np.exp(-xvalsBF) -1

    sdFig, ax = plt.subplots()
    ax.errorbar(zvals_sn, dL_obs/zvals_sn, yerr=error/zvals_sn, label=r"$\mathrm{Observation}$", fmt=".",  ecolor=Colors["d_L_obs"], capsize=4, elinewidth=1.5, color="red", markersize=7)
    ax.plot(zvals_pred, dL_pred/zvals_pred, color=Colors["d_L_fid"], label=r"$\mathrm{Fiducial\ cosmology}$")
    ax.plot(zvalsBF, dLBF/zvalsBF, color=Colors["d_L_best"], label=r"$\mathrm{Best\ fit}$")
    ax.set_xscale("log")
    ax.set_ylim(3.5,8)
    ax.set_xlim(0.005,1.45)
    ax.set_xlabel(lbls["z"])
    ax.set_ylabel(r"$\mathrm{Gpc}$")
    ax.set_title(r"$\mathrm{Luminosity\ distance},\ d_L$", loc="left")
    ax.legend(loc="upper left", fancybox=True)  
    ax.minorticks_on()

    
    save_push(sdFig, "supernova_data")


def goodness_of_fit():
    N = len(Sdata["d_L"])
    chi2 = Sfit["chi2"]
    chi2_over_N = chi2/N

    
    counts, bins = np.histogram(chi2_over_N, bins=150)
    sigma, mu = np.std(chi2_over_N), np.mean(chi2_over_N)

    gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))

    # goodnes, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(12,12))
    goodnes, ax1 = plt.subplots()
    hstg = ax1.hist(bins[:-1], bins, weights=counts, color=Colors["hist"], label=r"$\mathrm{samples}$", density=True)
    ax1.plot(bins, gaussian, color=Colors["gaussian"], label=r"$\mathrm{fit}$")
    ax1.axvline(mu, ls="--", color="black", lw=2, label=r"$\mu={mu:.3f}$".format(mu=mu))
    ax1.fill_between(bins, 0, gaussian, where=sigma > abs(bins-mu), color="blue", alpha=0.2, label=r"$\mu\pm 1\sigma$; $\sigma={sigma:.3f}$".format(sigma=sigma))
    # ax1.text(1.4, 2, r"$\mu = {mu:.3f} \\ \sigma = {sigma:.3f}$".format(mu=mu, sigma=sigma))
    ax1.set_title(r"$\mathrm{Goodness\ of\ fit}$", loc="left")
    ax1.set_xlabel(r"$\chi^2/N$")
    ax1.legend(loc="best", fancybox=True)



    # oneSDthres = 3.53
    # chi2min = np.min(chi2)
    # accepted = (chi2 - chi2min) < oneSDthres
    # chi2 = np.where(accepted, chi2, np.nan)
    # chi2 = chi2[np.isfinite(chi2)]
    # chi2_over_N = chi2/N
    
    # counts, bins = np.histogram(chi2_over_N, bins=150)
    # sigma, mu = np.std(chi2_over_N), np.mean(chi2_over_N)

    # gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))
    # hstg = ax2.hist(bins[:-1], bins, weights=counts, color=Colors["hist"], label="Samples", density=True)
    # ax2.plot(bins, gaussian, color=Colors["gaussian"], label="Fitted pdf")
    # ax2.axvline(mu, ls="--", color="black", lw=2)
    # ax2.fill_between(bins, 0, gaussian, where=sigma > abs(bins-mu), color="blue", alpha=0.2)


    save_push(goodnes, "goodnes_of_fit")





def omega_restrictions_plot():
    """Make a plot of the 1SD confidence of the chi2 values in the OmegaM, OmegaLambda plane.
    """
    chi2 = Sfit["chi2"]
    OmegaM = Sfit["OmegaM"]
    OmegaK = Sfit["OmegaK"]
    # OmegaLambda = 1-(OmegaK+OmegaM)


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
    ax.set_title(r"$1\sigma\ \mathrm{confidence\ plot\ of}\ (\Omega_M \;\mathrm{x}\; \Omega_\Lambda)$", loc="left")
    ax.minorticks_on()

    # fig.savefig(plot_path+"omega_restrictions.pdf", bbox_inches=None)

    save_push(OmegaPlane, "omega_plane")

def prob_plots():
    OmegaM = Sfit["OmegaM"]
    OmegaK = Sfit["OmegaK"]

    countM, binsM = np.histogram(OmegaM, bins=150)
    countK, binsK = np.histogram(OmegaK, bins=150)


    sigmaM, muM = np.std(OmegaM), np.mean(OmegaM)
    sigmaK, muK = np.std(OmegaK), np.mean(OmegaK)

    gaussianM = 1/(sigmaM*np.sqrt(2*np.pi))*np.exp(-(binsM-muM)**2/(2*sigmaM**2))
    gaussianK = 1/(sigmaK*np.sqrt(2*np.pi))*np.exp(-(binsK-muK)**2/(2*sigmaK**2))

    probs, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(12,12))
    histM = ax1.hist(binsM[:-1], binsM, weights=countM, color=Colors["hist"], label=r"$\mathrm{samples}$", density=True)
    gausM, = ax1.plot(binsM, gaussianM, color=Colors["gaussian"], label=r"$\mathrm{fit}$")

    histK = ax2.hist(binsK[:-1], binsK, weights=countK, color=Colors["hist"], density=True)
    gausK, = ax2.plot(binsK, gaussianK, color=Colors["gaussian"])

    muMl = ax1.axvline(muM, ls="--", color="black", lw=2, label=r"$\mu$")
    mumK = ax2.axvline(muK, ls="--", color="black", lw=2)

    oneSM = ax1.fill_between(binsM, 0, gaussianM, where=sigmaM > abs(binsM-muM), color="blue", alpha=0.2, label=r"$\mu\pm 1\sigma$")
    oneSK = ax2.fill_between(binsK, 0, gaussianK, where=sigmaK > abs(binsK-muK), color="blue", alpha=0.2)

    ax1.set_xlabel(lbls["OmegaM"])
    ax2.set_xlabel(lbls["OmegaK"])
    # ax1.set_ylabel("Probability")
    # ax2.set_ylabel("Probability")
    ax1.text(0.45, 3, r"$\mu_M = {mu:.3f} \\ \sigma_M = {sigma:.3f}$".format(mu=muM, sigma=sigmaM))
    ax2.text(-.9, 1, r"$\mu_k = {mu:.3f} \\ \sigma_k = {sigma:.3f}$".format(mu=muK, sigma=sigmaK))
    probs.suptitle(r"$\mathrm{Posterior\ pdf\ of}\ \Omega_M\ \mathrm{and}\ \Omega_k$", x=.1, horizontalalignment="left")
    probs.legend(loc="center right", fancybox=True, bbox_to_anchor=[0.97, 0.90], ncol=4, fontsize=26)
    ax1.minorticks_on()
    ax2.minorticks_on()



    save_push(probs, "probs_M_K.pdf")
    print(f"For OmegaM:")
    print(f"mu = {muM:.3f}, sigma = {sigmaM:.3f}")
    print(f"For OmegaK:")
    print(f"mu = {muK:.3f}, sigma = {sigmaK:.3f}")


def posterior_pdf():
    h = Sfit["h"]

    #   Make H0 from h
    # H0 =100 * h
    H0 = h

    #create count and bins
    counts, bins = np.histogram(H0, bins=150)

    # Gaussian func
    sigma, mu = np.std(H0), np.mean(H0)
    gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))

    ppdf, ax = plt.subplots()
    hstg = ax.hist(bins[:-1], bins, weights=counts, color=Colors["hist"], label="Samples", density=True)
    ax.plot(bins, gaussian, color=Colors["gaussian"], label="Fitted pdf")
    ax.axvline(mu, ls="--", color="black", lw=2)
    ax.fill_between(bins, 0, gaussian, where=sigma > abs(bins-mu), color="blue", alpha=0.2)
    ax.text(.71, 50, r"$\mu = {mu:.3f} \\ \sigma = {sigma:.3f}$".format(mu=mu, sigma=sigma))
    ax.set_xlabel(r"$100\ \mathrm{kms}^{-1}\mathrm{Mpc}^{-1}$")
    # ax.set_ylabel("Probability")
    ax.set_title(r"$\mathrm{Posterior\ pdf\ of}\ H_0$", loc="left")
    ax.minorticks_on()

    # l1 = ax.legend(loc="best", fancybox=True)

    save_push(ppdf, "posterior_pdf")
    print(f"H0:")
    print(f"mu = {mu:.3f}, sigma = {sigma:.3f}")

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
    prob_plots()
    goodness_of_fit()
    omega_restrictions_plot()
    posterior_pdf()
    # create_table()

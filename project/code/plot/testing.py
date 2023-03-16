from plot_utils import *

Cosmology = Data("backgroundcosmology.csv")
ValueTab = Data("table_of_values.csv")

x_min = Cosmology["x"][0]
# embed()
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




if __name__=="__main__":
    testing_Omegas()
    testing_Hp()
    testing_eta()
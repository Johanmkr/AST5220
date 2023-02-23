from plot_utils import *

Cosmology = Data("backgroundcosmology.txt")

def testing_Omegas():
    """
    Plot for testing the Omegas evolution
    """
    xvals = Cosmology["x"]
    OmegaRad = Cosmology["OmegaR"] + Cosmology["OmegaNu"]
    OmegaM = Cosmology["OmegaB"] + Cosmology["OmegaCDM"]
    
    omegaFig, ax1 = plt.subplots()
    ax1.plot(xvals, Cosmology["OmegaLambda"]+OmegaRad+OmegaM, label="Sum", color="black", ls="--")
    ax1.plot(xvals, OmegaRad, label=r"$\Omega_\mathrm{rad}$", color="orange")
    ax1.plot(xvals, OmegaM, label=r"$\Omega_\mathrm{M}$", color="green")
    ax1.plot(xvals, Cosmology["OmegaLambda"], label=r"$\Omega_\Lambda$", color="purple")
    ax1.set_xlabel(r"$x$")
    ax1.set_title(r"Density fractions $\Omega_X$", loc="left")
    omegaFig.legend(loc='center right', bbox_to_anchor=(.95, 0.5),
          ncol=1, fancybox=True)
    
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
    line1, = ax1.plot(xvals, ratioddHp, color="red", label=r"$\frac{1}{\mathcal{H}}\frac{\mathrm{d}^2\mathcal{H}}{\mathrm{d}x^2}$")
    line2, = ax1.plot(xvals, ratiodHp, color="blue", label=r"$\frac{1}{\mathcal{H}}\frac{\mathrm{d}\mathcal{H}}{\mathrm{d}x}$")
    rad_area = ax1.axvspan(-20, -11, color="orange", alpha=.15, label=r"$\Omega_\mathrm{rad}$")
    mat_area = ax1.axvspan(-5, -1.5, color="green", alpha=.15, label=r"$\Omega_\mathrm{M}$")
    lam_area = ax1.axvspan(1.5, 5, color="purple", alpha=.15, label=r"$\Omega_\Lambda$")
    ax1.set_title(r"Sanity check of $\mathcal{H}(x)$", loc="left")
    ax1.set_xlabel(r"$x$")
    legend1 = ax1.legend([line1, line2], [line1.get_label(), line2.get_label()], loc="center left", fancybox=True)
    legend2 = HpFig.legend([rad_area, mat_area, lam_area], [rad_area.get_label(), mat_area.get_label(), lam_area.get_label()], loc="upper right", fancybox=True, ncols=3, bbox_to_anchor=[0.99,1])

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
    ax1.plot(xvals, etaHpc, color="blue", label=r"$\frac{\eta\mathcal{H}}{c}$")
    ax1.plot(xvals, detaHpc, color="red", label=r"$\frac{\mathcal{H}}{c}\frac{\mathrm{d}\eta}{\mathrm{d}x}$")
    ax1.set_title(r"Sanity check for $\eta(x)$", loc="left")
    ax1.set_xlabel(r"$x$")
    # ax1.set_yscale("log")
    ax1.legend(loc="best", fancybox=True)

    save_push(etaFig, "eta_test")




if __name__=="__main__":
    testing_Omegas()
    testing_Hp()
    testing_eta()
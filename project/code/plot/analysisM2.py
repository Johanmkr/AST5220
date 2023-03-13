from plot_utils import *

Recombination = Data("recombination.csv")
# embed()



def plot_Xe():
    xvals = Recombination["x"]
    Xe_vals = Recombination["Xe"]

    Xe, ax = plt.subplots()
    ax.plot(xvals, Xe_vals)
    ax.set_yscale("log")
    # ax.set_xlim(-12,0)

    save_push(Xe, "Xe_electron_fraction")



def plot_tau():
    xvals = Recombination["x"]
    tau_vals = Recombination["tau"]
    dtauvals = Recombination["dtaudx"]
    ddtauvals = Recombination["ddtauddx"]


    Tau, ax = plt.subplots()
    ax.plot(xvals, tau_vals, color="blue", label="tau")
    ax.plot(xvals, -dtauvals, color="red", label="dtau")
    ax.plot(xvals, ddtauvals, color="green", label="ddtau")
    ax.set_yscale("log")

    ax.legend(loc="best")

    save_push(Tau, "tau_of_x")

def plot_g():
    xvals = Recombination["x"]
    g_vals = Recombination["g"]
    dgvals = Recombination["dgdx"]
    ddgvals = Recombination["ddgddx"]
    # embed()

    G1, ax1 = plt.subplots()
    ax1.plot(xvals, g_vals, color="blue", label="g")
    ax1.legend()
    save_push(G1, "g_of_x")

    G2, ax2 = plt.subplots()
    ax2.plot(xvals, dgvals, color="blue", label="dgdx")
    ax2.legend()
    save_push(G2, "dg_of_x")

    G3, ax3 = plt.subplots()
    ax3.plot(xvals, ddgvals, color="blue", label="ddgddx")
    ax3.legend()
    save_push(G3, "ddg_of_x")


if __name__=="__main__":
    plot_Xe()
    plot_tau()
    plot_g()
from plot_utils import *

class Recombination(Milestone):
    def __init__(self, data:Data)->None:
        super().__init__(data)

    def plot_Xe(self)->None:
        XeFig, ax = plt.subplots()

        # Plot Xe from Saha and peebles
        XeLine, = ax.plot(self.x, self.Xe, color=Colors["Xe"])
        XeLine.set_label("Peebles")

        # Plot Xe from Saha only
        SahaLine, = ax.plot(self.x, self.XeSaha, color=Colors["XeSaha"], ls="dashed", lw=2)
        SahaLine.set_label("Saha")

        # misc
        ax.set_xlim(-8,-4)
        ax.set_ylim(0.0001, 1.2)
        ax.set_yscale("log")
        ax.set_xlabel(lbls["x"])
        ax.set_ylabel(lbls["Xe"])
        ax.set_title(r"Free electron fraction, $X_e$", loc="left")

        # Make legend
        L1 = ax.legend(loc="best", fancybox=True)
        
        # Make self variable and create plot
        self.XeFig = XeFig
        save_push(self.XeFig, "Xe_plot")

    def plot_tau(self)->None:
        fig, ax = plt.subplots()

        # Make lines
        tau, = ax.plot(self.x, self.tau, label=lbls["tau"], color=Colors["tau"])
        dtau, = ax.plot(self.x, -self.dtaudx, label=lbls["-dtaudx"], color=Colors["dtaudx"], ls="dashed", lw=2)
        ddtau, = ax.plot(self.x, self.ddtauddx, label=lbls["ddtauddx"], color=Colors["ddtauddx"], ls="dashdot", lw=2)

        # Titles
        ax.set_title(r"Optical depth, $\tau$", loc="left")
        ax.set_xlabel(lbls["x"])


        # Misc
        ax.set_yscale("log")
        ax.set_xlim(-10,-4)
        ax.set_ylim(1e-5, 1e5)

        # Make legend
        L1 = ax.legend(loc="best", fancybox=True)

        # Make self variable and create plot
        self.TauFig = fig
        save_push(self.TauFig, "optical_depth")

    def plot_g(self)->None:
        fig, ax = plt.subplots()

        # Make lines
        g, = ax.plot(self.x, self.g/np.abs(self.g).max(), label=lbls["g"], color=Colors["g"])
        dg, = ax.plot(self.x, self.dgdx/np.abs(self.dgdx).max(), label=lbls["dgdx"], color=Colors["dgdx"], ls="dashed", lw=2)
        ddg, = ax.plot(self.x, self.ddgddx/np.abs(self.ddgddx).max(), label=lbls["ddgddx"], color=Colors["ddgddx"], ls="dashdot", lw=2)

        # Titles
        ax.set_xlabel(lbls["x"])
        ax.set_title(r"Visibility function, $\tilde{g}$", loc="left")

        # Misc
        ax.set_xlim(-7.5,-6)

        # Make legend
        L1 = ax.legend(loc="best", fancybox=True)

        # Make self variable and create plot
        self.GFig = fig
        save_push(self.GFig, "visibility_function")


    def make_plots(self):
        self.plot_Xe()
        self.plot_tau()
        self.plot_g()
        





# def plot_Xe():
#     xvals = Recombination["x"]
#     Xe_vals = Recombination["Xe"]

#     Xe, ax = plt.subplots()
#     ax.plot(xvals, Xe_vals)
#     ax.set_yscale("log")
#     ax.set_xlim(-12,0)

#     save_push(Xe, "Xe_electron_fraction")



# def plot_tau():
#     xvals = Recombination["x"]
#     tau_vals = Recombination["tau"]
#     dtauvals = Recombination["dtaudx"]
#     ddtauvals = Recombination["ddtauddx"]


#     Tau, ax = plt.subplots()
#     ax.plot(xvals, tau_vals, color="blue", label="tau")
#     ax.plot(xvals, -dtauvals, color="red", label="dtau")
#     ax.plot(xvals, ddtauvals, color="green", label="ddtau")
#     ax.set_yscale("log")
#     ax.set_xlim(-12,0)


#     ax.legend(loc="best")

#     save_push(Tau, "tau_of_x")

# def plot_g():
#     xvals = Recombination["x"]
#     g_vals = Recombination["g"]
#     dgvals = Recombination["dgdx"]
#     ddgvals = Recombination["ddgddx"]
#     # embed()

#     G1, ax1 = plt.subplots()
#     ax1.plot(xvals, g_vals, color="blue", label="g")
#     ax1.set_xlim(-12,0)
#     ax1.legend()
#     save_push(G1, "g_of_x")

#     G2, ax2 = plt.subplots()
#     ax2.plot(xvals, dgvals, color="blue", label="dgdx")
#     ax2.set_xlim(-12,0)
#     ax2.legend()
#     save_push(G2, "dg_of_x")

#     G3, ax3 = plt.subplots()
#     ax3.plot(xvals, ddgvals, color="blue", label="ddgddx")
#     ax3.set_xlim(-12,0)
#     ax3.legend()
#     save_push(G3, "ddg_of_x")


if __name__=="__main__":
    # plot_Xe()
    # plot_tau()
    # plot_g()
    Rec = Recombination(Data("recombination.csv"))
    Rec.make_plots()
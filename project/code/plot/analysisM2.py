from plot_utils import *

class Recombination(Milestone):
    def __init__(self, data:Data)->None:
        super().__init__(data)
        self.AnalData = Data("recomb_analysis.csv")
        self.x_LS, self.x_rec, self.x_recSaha = self.AnalData["x"][0], self.AnalData["x"][1], self.AnalData["x"][2]
        # self.t_LS, self.t_rec, self.t_recSaha = self.AnalData["t"][0], self.AnalData["t"][1], self.AnalData["t"][2]
        # self.r_s_LS, self.r_s_rec, self.r_s_recSaha = self.AnalData["r_s"][0], self.AnalData["r_s"][1], self.AnalData["r_s"][2]

    def plot_Xe(self)->None:
        XeFig, ax = plt.subplots()

        # Plot Xe from Saha and peebles
        XeLine, = ax.plot(self.x, self.Xe, color=Colors["Xe"])
        XeLine.set_label(lbls["Xe"] + " (Peebles)")

        # Plot Xe from Saha only
        SahaLine, = ax.plot(self.x, self.XeSaha, color=Colors["XeSaha"], ls="dashed", lw=2)
        SahaLine.set_label(lbls["Xe"] + " (Saha)")

        # Misc
        ax.set_xlim(-8,-4)
        ax.set_ylim(0.0001, 1.2)
        ax.set_yscale("log")
        ax.set_xlabel(lbls["x"])
        # ax.set_ylabel(lbls["Xe"])
        ax.set_title(r"Free electron fraction, $X_e$", loc="left")

        # Plot lines
        ax.axvline(self.x_rec, lw=2, ls="dashdot", color="black", label="Recombination")
        ax.axvline(self.x_recSaha, lw=2, ls="dotted", color="black", label="Saha recombination")
        ax.axhline(self.Xe[np.argmin(np.abs(self.x))], lw=2, ls="dashed", color="chocolate", label="Freeze out")

        # Make legend
        ax.minorticks_on()
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
        lines = [tau, dtau, ddtau]

        # Titles
        ax.set_title(r"Optical depth, $\tau$", loc="left")
        ax.set_xlabel(lbls["x"])


        # Misc
        ax.set_yscale("log")
        ax.set_xlim(-10,-4)
        ax.set_ylim(1e-5, 1e5)
        lsline = ax.axvline(self.x_LS, lw=2, ls="dashed", color="black", label="Last scattering")

        # Make legend
        ax.minorticks_on()
        L1 = ax.legend(lines, [line.get_label() for line in lines], loc="best", fancybox=True)
        L2 = fig.legend([lsline], ["Last Scattering"], loc="upper right", fancybox=True)

        # Make self variable and create plot
        self.TauFig = fig
        save_push(self.TauFig, "optical_depth")

    def plot_g(self)->None:
        fig, ax = plt.subplots()

        # Make lines
        g, = ax.plot(self.x, self.g, label=lbls["g"], color=Colors["g"])
        dg, = ax.plot(self.x, self.dgdx/10, label=lbls["dgdx10"], color=Colors["dgdx"], ls="dashed", lw=2)
        ddg, = ax.plot(self.x, self.ddgddx/300, label=lbls["ddgddx300"], color=Colors["ddgddx"], ls="dashdot", lw=2)
        lines = [g, dg, ddg]

        # Titles
        ax.set_xlabel(lbls["x"])
        ax.set_title(r"Visibility function, $\tilde{g}$", loc="left")

        # Misc
        ax.set_xlim(-7.5,-6)
        lsline = ax.axvline(self.x_LS, lw=2, ls="dashed", color="black")


        # Make legend
        ax.minorticks_on()
        L1 = ax.legend(lines, [line.get_label() for line in lines], loc="upper right", fancybox=True)
        L2 = fig.legend([lsline], ["Last Scattering"], loc="upper right", fancybox=True)

        # Make self variable and create plot
        self.GFig = fig
        save_push(self.GFig, "visibility_function")

    def make_table(self)->None:
        styler = self.AnalData.dF.style
        styler.format({"x": '{:.4f}', "z": '{:.2f}', "t [Myr]": '{:.4f}', "r_s [Mpc]": '{:.2f}'})
        styler.hide(axis="index")
        styler.hide_columns([self.AnalData.dF.keys()[-1]])
        styler.to_latex(latex_path + "recomb_analysis.tex")


    def make_plots(self):
        self.plot_Xe()
        self.plot_tau()
        self.plot_g()
        




if __name__=="__main__":
    Rec = Recombination(Data("recombination.csv"))
    Rec.make_plots()
    # Rec.make_table()
    # print(Rec.AnalData)
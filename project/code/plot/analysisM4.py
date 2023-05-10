from plot_utils import *


class PowerSpectrum:
    def __init__(self):
        # Make data objects
        self.Cell = Milestone(Data("cellss.csv"))
        self.Theta = Milestone(Data("theta_l.csv"))
        self.Cell_sep = Milestone(Data("cell_separated.csv"))
        self.bessel = Milestone(Data("bessel.csv"))
        self.MPS = Milestone(Data("mps.csv"))
        self.CL_integrand = Milestone(Data("cl_integrand.csv"))
        self.LOS_integrand = Milestone(Data("LOS_integrand.csv"))
        self.figures = {}
    
    # def add_figure(self, figname, **kwargs):
    #     fig, ax = plt.subplots(**kwargs)
    #     self.figures[figname] = fig
    #     return fig
        
    def transfer_function(self):
        info_dicts = []
        ell_values = [6, 100, 200, 500, 1000]
        T_values = []

        #   Make information about each line and figure
        for i in range(len(ell_values)):
            info_dict = dict(
                lw = THINLINE,
                label = r"$l={ell}$".format(ell=ell_values[i])
            )
            info_dicts.append(info_dict)
            T_values.append(f"self.Theta.T_{ell_values[i]}")

        ax_setter_info = dict(
            xlabel=r"$c/H_0$",
            title=r"$\mathrm{Transfer\ function\ }\Theta_l$",
            xlim=[0,500],
            ylim=[-0.005, 0.015]
        )

        ThetaPlot = MAKEPLOT("transfer_function")
        for i in range(len(ell_values)):
            ThetaPlot.plot_line(self.Theta.ckH, eval(T_values[i]), **info_dicts[i])
        
        ThetaPlot.set_ax_info(**ax_setter_info)
        ThetaPlot.set_minor_ticks()
        ThetaPlot.set_legend(fancybox=True, loc="upper right", ncols=2)

        self.figures["ThetaPlot"] = ThetaPlot

        ThetaPlot.finished()

    def C_l_integrand(self):
        info_dicts = []
        ell_values = [6, 100, 200, 500, 1000]
        T2_values = []

        #   Make information about each line and figure
        for i in range(len(ell_values)):
            info_dict = dict(
                lw = THINLINE,
                label = r"$l={ell}$".format(ell=ell_values[i])
            )
            info_dicts.append(info_dict)
            T2_values.append(f"self.CL_integrand.T2_{ell_values[i]}")

        Cl_int_plot = MAKEPLOT("C_l_integrand")
        for i in range(len(ell_values)):
            Cl_int_plot.plot_line(self.CL_integrand.k, eval(T2_values[i]), **info_dicts[i])

        ax_setter_info = dict(
            xlabel=r"$c/H_0$",
            title=r"$\mathrm{Integrand\ of\ Power\ Spectrum:\ }|\Theta_l(k)|^2/k$",
            ylabel=r"$H_0/c$",
            xlim=[0,500],
            ylim=[0,5e-7]
        )

        Cl_int_plot.set_ax_info(**ax_setter_info)
        Cl_int_plot.set_minor_ticks()
        Cl_int_plot.set_legend(fancybox=True, loc="upper right", ncols=2)

        self.figures["Cl_int_plot"] = Cl_int_plot

        Cl_int_plot.finished()

    def PowerSpectrum_plot(self):
        subparts = ["SW", "ISW", "DOP", "POL"]
        info_dicts = []
        sub_evals = []

        for i in range(len(subparts)):
            info_dict = dict(
                lw=THINLINE,
                label = r"$\mathrm{%s}$" % subparts[i],
                ls="dashed"
            )
            info_dicts.append(info_dict)
            sub_evals.append(f"self.Cell_sep.cell_{subparts[i]}")
        
        ax_setter_info = dict(
            xlabel = r"$l$",
            ylabel = r"$\frac{l(l+1)}{2\pi}(10^6T_\mathrm{CMB0})^2$",
            title = r"$\mathrm{Power\ spectrum\ for\ radiation,\ } C_l$",
            xscale = "log"
        )

        PSplot = MAKEPLOT("power_spectrum")
        PSplot.plot_line(self.Cell.ell, self.Cell.cell_TT, label=r"$C_l$", color="blue")

        for i in range(len(subparts)):
            PSplot.plot_line(self.Cell_sep.ell, eval(sub_evals[i]), **info_dicts[i])
        
        PSplot.set_ax_info(**ax_setter_info)
        PSplot.set_minor_ticks()
        PSplot.set_legend(fancybox=True, loc="best", ncols=2)

        self.figures["PSplot"] = PSplot

        PSplot.finished()

    def MPS_plot(self):


        ax_setter_info = dict(
            xscale="log",
            yscale="log",
            title=r"$\mathrm{Matter\ power\ spectrum,\ } P(k)$",
            xlabel=r"$h/\mathrm{Mpc}$",
            ylabel=r"$(\mathrm{Mpc}/h)^3$"
        )
        MPSplot = MAKEPLOT("matter_power_spectrum")
        MPSplot.plot_line(self.MPS.k, self.MPS.Pk, label=r"$P(k)$", color="blue")
        MPSplot.set_ax_info(**ax_setter_info)

        MPSplot.finished()



        
    


# power = Milestone(Data("cellss.csv"))




# power = Milestone(Data("cellss.csv"))
# Theta = Milestone(Data("theta_l.csv"))


# fig, ax = plt.subplots()
# LW=1.
# ax.plot(Theta.ckH, Theta.T_6, lw=LW, label="l=6")
# ax.plot(Theta.ckH, Theta.T_100, lw=LW, label="l=100")
# ax.plot(Theta.ckH, Theta.T_200, lw=LW, label="l=200")
# ax.plot(Theta.ckH, Theta.T_500, lw=LW, label="l=500")
# ax.plot(Theta.ckH, Theta.T_1000, lw=LW, label="l=1000")
# ax.set_xlim(0,500)
# ax.set_ylim(-0.01, 0.01)
# ax.legend()


# fig, ax = plt.subplots()
# ax.plot(power.ell, power.cell_TT)
# ax.set_xscale("log")


# plt.show()

if __name__=="__main__":
    PS = PowerSpectrum()
    PS.transfer_function()
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
    
    def add_figure(self, figname, **kwargs):
        fig, ax = plt.subplots(**kwargs)
        self.figures[figname] = fig
        return fig
        
    def transfer_function(self):
        

        
    


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
from plot_utils import *


# class PowerSpectrum:
#     def __init__(self)
#         self.Cell = Milestone(Data("cellss.csv"))
#         self.Theta = Milestone(Data("theta_l.csv"))
#         self.figures = []
    
    


# power = Milestone(Data("cellss.csv"))

power = Milestone(Data("cellss.csv"))
Theta = Milestone(Data("theta_l.csv"))


fig, ax = plt.subplots()
LW=1.
ax.plot(Theta.ckH, Theta.T_6, lw=LW, label="l=6")
ax.plot(Theta.ckH, Theta.T_100, lw=LW, label="l=100")
ax.plot(Theta.ckH, Theta.T_200, lw=LW, label="l=200")
ax.plot(Theta.ckH, Theta.T_500, lw=LW, label="l=500")
ax.plot(Theta.ckH, Theta.T_1000, lw=LW, label="l=1000")
ax.set_xlim(0,500)
ax.set_ylim(-0.01, 0.01)
ax.legend()


fig, ax = plt.subplots()
ax.plot(power.ell, power.cell_TT)
ax.set_xscale("log")


plt.show()
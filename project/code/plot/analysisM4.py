from plot_utils import *
import matplotlib as mpl
mpl.use("Qt5Agg")


power = Milestone(Data("cellss.csv"))


fig, ax = plt.subplots()

ax.plot(power.ell, power.cell_TT)

plt.show()
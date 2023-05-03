from plot_utils import *


power = Milestone(Data("cellss.csv"))


fig, ax = plt.subplots()

ax.plot(power.ell, power.cell_TT)

plt.show()
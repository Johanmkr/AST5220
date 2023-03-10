from plot_utils import *

Recombination = Data("recombination.csv")
# embed()



def plot_Xe():
    xvals = Recombination["x"]
    Xe_vals = Recombination["Xe"]

    Xe, ax = plt.subplots()
    ax.plot(xvals, Xe_vals)
    ax.set_yscale("log")

    save_push(Xe, "Xe_electron_fraction")



if __name__=="__main__":
    plot_Xe()
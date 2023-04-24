from plot_utils import *

# class Perturbation:
#     def __init__(self, list:DataList)->None:
#         self.DataList = DataList

class Perturbation:
    def __init__(self, files:list)->None:
        self.k1Data = Data(files[0])
        self.k2Data = Data(files[1])
        self.k3Data = Data(files[2])
        for key in self.k1Data.get_keys():
            exec(f"self.k1_{key} = self.k1Data[key]")
        for key in self.k2Data.get_keys():
            exec(f"self.k2_{key} = self.k2Data[key]")
        for key in self.k3Data.get_keys():
            exec(f"self.k3_{key} = self.k3Data[key]")




    def potentials_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_Phi, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, self.k2_Phi, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, self.k3_Phi, color=Colors["k3"], label=lbls["k3"])
        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "potentials")

    def monopole_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_T0, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, self.k2_T0, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, self.k3_T0, color=Colors["k3"], label=lbls["k3"])
        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "monopole")

    def dipole_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_T1, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, self.k2_T1, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, self.k3_T1, color=Colors["k3"], label=lbls["k3"])
        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "dipole")

    def quadrapole_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_T2, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, self.k2_T2, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, self.k3_T2, color=Colors["k3"], label=lbls["k3"])
        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "quadrapole")

    def delta_plot(self)->None:
        fig, ax = plt.subplots()
        # self.k1_delta_gamma = 4*self.k1_T0
        # self.k2_delta_gamma = 4*self.k2_T0
        # self.k3_delta_gamma = 4*self.k3_T0
        ax.plot(self.k1_x, self.k1_delta_cdm, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, self.k2_delta_cdm, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, self.k3_delta_cdm, color=Colors["k3"], label=lbls["k3"])
        ax.plot(self.k1_x, abs(self.k1_delta_b), color=Colors["k1"], ls="dashed")
        ax.plot(self.k2_x, abs(self.k2_delta_b), color=Colors["k2"], ls="dashed")
        ax.plot(self.k3_x, abs(self.k3_delta_b), color=Colors["k3"], ls="dashed")

        ax.set_yscale("log")
        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "delta")


    def velocity_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_v_cdm, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, self.k2_v_cdm, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, self.k3_v_cdm, color=Colors["k3"], label=lbls["k3"])
        ax.plot(self.k1_x, abs(self.k1_v_b), color=Colors["k1"], ls="dashed")
        ax.plot(self.k2_x, abs(self.k2_v_b), color=Colors["k2"], ls="dashed")
        ax.plot(self.k3_x, abs(self.k3_v_b), color=Colors["k3"], ls="dashed")

        ax.set_title(r"")

        ax.set_yscale("log")
        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "velocity")


    def make_plots(self)->None:
        # self.potentials_plot()
        # self.monopole_plot()
        # self.dipole_plot()
        self.delta_plot()
        self.velocity_plot()


if __name__=="__main__":
    Pert = Perturbation(["perturbations_k0.001.csv", "perturbations_k0.01.csv", "perturbations_k0.1.csv"])
    Pert.make_plots()
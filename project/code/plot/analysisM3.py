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



    def phi_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_Phi, label="k1")
        ax.plot(self.k2_x, self.k2_Phi, label="k2")
        ax.plot(self.k3_x, self.k3_Phi, label="k3")
        ax.minorticks_on()
        ax.legend()
        plt.show()

    def monopole_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_T0, label="k1")
        ax.plot(self.k2_x, self.k2_T0, label="k2")
        ax.plot(self.k3_x, self.k3_T0, label="k3")
        ax.minorticks_on()
        ax.legend()
        plt.show()

    def dipole_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_T1, label="k1")
        ax.plot(self.k2_x, self.k2_T1, label="k2")
        ax.plot(self.k3_x, self.k3_T1, label="k3")
        ax.minorticks_on()
        ax.legend()
        plt.show()

    def delta_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_delta_b, label="k1")
        ax.plot(self.k2_x, self.k2_delta_b, label="k2")
        ax.plot(self.k3_x, self.k3_delta_b, label="k3")
        ax.minorticks_on()
        ax.legend()
        plt.show()

    def velocity_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_v_b, label="k1")
        ax.plot(self.k2_x, self.k2_v_b, label="k2")
        ax.plot(self.k3_x, self.k3_v_b, label="k3")
        ax.minorticks_on()
        ax.legend()
        plt.show()

    def make_plots(self)->None:
        self.phi_plot()
        self.monopole_plot()
        self.dipole_plot()
        self.delta_plot()
        self.velocity_plot()


if __name__=="__main__":
    Pert = Perturbation(["perturbations_k0.001.csv", "perturbations_k0.01.csv", "perturbations_k0.1.csv"])
    Pert.make_plots()
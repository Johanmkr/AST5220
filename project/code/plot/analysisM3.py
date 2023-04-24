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
        self.LWpoles = 2

    def set_k_lables(self, fig, lines):
        leg = fig.legend(lines, [line.get_label() for line in lines], loc="upper right", fancybox=True, ncol=3, title=r"$k\ [\mathrm{Mpc}^{-1}]=$", fontsize=20, title_fontsize=20, bbox_to_anchor=[0.9675, 1.015])
        leg._legend_box.align = "left"


    def potentials_plot(self)->None:
        fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(12,12), sharex=True)
        ax1.plot(self.k1_x, self.k1_Phi, color=Colors["k1"], label=lbls["k1"])
        ax1.plot(self.k2_x, self.k2_Phi, color=Colors["k2"], label=lbls["k2"])
        ax1.plot(self.k3_x, self.k3_Phi, color=Colors["k3"], label=lbls["k3"])
        ax2.plot(self.k1_x, 1e4*(self.k1_Phi+self.k1_Psi), color=Colors["k1"])
        ax2.plot(self.k2_x, 1e4*(self.k2_Phi+self.k2_Psi), color=Colors["k2"])
        ax2.plot(self.k3_x, 1e4*(self.k3_Phi+self.k3_Psi), color=Colors["k3"])

        ax2.set_xlabel(lbls["x"])
        ax1.set_title(r"$\mathrm{Potential}\ \Phi$", loc="left")
        ax2.set_title(r"$\mathrm{Sum:\ }\ 10^4\cdot(\Phi+\Psi)$", loc="left")

        ax2.minorticks_on()
        ax1.legend(loc="best", fancybox=True)
        
        save_push(fig, "potentials")

    def monopole_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, 4*self.k1_T0, lw=self.LWpoles, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, 4*self.k2_T0, lw=self.LWpoles, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, 4*self.k3_T0, lw=self.LWpoles, color=Colors["k3"], label=lbls["k3"])

        ax.set_xlabel(lbls["x"])
        ax.set_title(r"$\mathrm{Overdensity}\ \delta_\gamma=4\Theta_0$", loc="left")
        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "monopole")

    def dipole_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, -3*self.k1_T1, lw=self.LWpoles, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, -3*self.k2_T1, lw=self.LWpoles, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, -3*self.k3_T1, lw=self.LWpoles, color=Colors["k3"], label=lbls["k3"])
        
        ax.set_xlabel(lbls["x"])
        ax.set_title(r"$\mathrm{Velocity}\ v_\gamma=-3\Theta_1$", loc="left")
        
        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "dipole")

    def quadrapole_plot(self)->None:
        fig, ax = plt.subplots()
        ax.plot(self.k1_x, self.k1_T2, lw=self.LWpoles, color=Colors["k1"], label=lbls["k1"])
        ax.plot(self.k2_x, self.k2_T2, lw=self.LWpoles, color=Colors["k2"], label=lbls["k2"])
        ax.plot(self.k3_x, self.k3_T2, lw=self.LWpoles, color=Colors["k3"], label=lbls["k3"])

        ax.set_xlabel(lbls["x"])
        ax.set_title(r"$\mathrm{Quadrapole}\ \Theta_2$", loc="left")

        ax.minorticks_on()
        ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "quadrapole")

    def delta_plot(self)->None:
        fig, ax = plt.subplots()
        line1, = ax.plot(self.k1_x, self.k1_delta_cdm, color=Colors["k1"], label=lbls["k1"])
        line2, = ax.plot(self.k2_x, self.k2_delta_cdm, color=Colors["k2"], label=lbls["k2"])
        line3, = ax.plot(self.k3_x, self.k3_delta_cdm, color=Colors["k3"], label=lbls["k3"])
        ax.plot(self.k1_x, abs(self.k1_delta_b), lw=2, color=Colors["k1"], ls="dashed")
        ax.plot(self.k2_x, abs(self.k2_delta_b), lw=2, color=Colors["k2"], ls="dashed")
        ax.plot(self.k3_x, abs(self.k3_delta_b), lw=2, color=Colors["k3"], ls="dashed")
        ax.set_title(r"$\mathrm{Overdensities\ } \delta_c\ \mathrm{and}\ \delta_b$", loc="left")

        line_lab, = ax.plot(np.nan, np.nan, color="grey")
        dash_lab, = ax.plot(np.nan, np.nan, color="grey", ls="dashed")
        ax.set_xlabel(lbls["x"])
        # ax.set_ylim(1e-1, 1e5)
        ax.set_yscale("log")
        ax.minorticks_on()
        
        lines=[line1, line2, line3]
        self.set_k_lables(fig, lines)
        ax.legend([line_lab, dash_lab], [r"$\delta_c$", r"$\delta_b$"], fancybox=True, loc="best")
        
        # ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "delta")


    def velocity_plot(self)->None:
        fig, ax = plt.subplots()
        line1, = ax.plot(self.k1_x, self.k1_v_cdm, color=Colors["k1"], label=lbls["k1"])
        line2, = ax.plot(self.k2_x, self.k2_v_cdm, color=Colors["k2"], label=lbls["k2"])
        line3, = ax.plot(self.k3_x, self.k3_v_cdm, color=Colors["k3"], label=lbls["k3"])
        ax.plot(self.k1_x, abs(self.k1_v_b), lw=2, color=Colors["k1"], ls="dashed")
        ax.plot(self.k2_x, abs(self.k2_v_b), lw=2, color=Colors["k2"], ls="dashed")
        ax.plot(self.k3_x, abs(self.k3_v_b), lw=2, color=Colors["k3"], ls="dashed")

        ax.set_title(r"$\mathrm{Velocities\ } v_c\ \mathrm{and}\ v_b$", loc="left")

        ax.set_xlabel(lbls["x"])
        ax.set_xlim(-17.5, 1)
        ax.set_yscale("log")
        ax.minorticks_on()
        # ax.legend(loc="best", fancybox=True)
        line_lab, = ax.plot(np.nan, np.nan, color="grey")
        dash_lab, = ax.plot(np.nan, np.nan, color="grey", ls="dashed")
        lines=[line1, line2, line3]
        self.set_k_lables(fig, lines)
        ax.legend([line_lab, dash_lab], [r"$v_c$", r"$v_b$"], fancybox=True, loc="best")
        
        save_push(fig, "velocity")


    def make_plots(self)->None:
        # self.potentials_plot()
        # self.monopole_plot()
        # self.dipole_plot()
        # self.quadrapole_plot()
        self.delta_plot()
        self.velocity_plot()


if __name__=="__main__":
    Pert = Perturbation(["perturbations_k0.001.csv", "perturbations_k0.01.csv", "perturbations_k0.1.csv"])
    Pert.make_plots()
from plot_utils import *

# class Perturbation:
#     def __init__(self, list:DataList)->None:
#         self.DataList = DataList
RECTOL = 0.15

class Perturbation:
    def __init__(self, files:list, x_rec=None, x_RM=None)->None:
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
        self.MarkLines = 2
        self.x_rec = x_rec
        self.x_RM = x_RM

    def set_k_lables(self, fig, lines):
        leg = fig.legend(lines, [line.get_label() for line in lines], loc="upper right", fancybox=True, ncol=3, title=r"$k\ [\mathrm{Mpc}^{-1}]=$", fontsize=20, title_fontsize=20, bbox_to_anchor=[0.95, 1.01])
        leg._legend_box.align = "left"

    def potentials_plot(self)->None:
        fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(12,12), sharex=True)
        line1, = ax1.plot(self.k1_x, self.k1_Phi, color=Colors["k1"], label=lbls["k1"])
        line2, = ax1.plot(self.k2_x, self.k2_Phi, color=Colors["k2"], label=lbls["k2"])
        line3, = ax1.plot(self.k3_x, self.k3_Phi, color=Colors["k3"], label=lbls["k3"])
        ax2.plot(self.k1_x, (self.k1_Phi+self.k1_Psi), color=Colors["k1"])
        ax2.plot(self.k2_x, (self.k2_Phi+self.k2_Psi), color=Colors["k2"])
        ax2.plot(self.k3_x, (self.k3_Phi+self.k3_Psi), color=Colors["k3"])

        ax1.set_xlim(-15,0)
        ax2.set_xlim(-15,0)
        ax2.set_ylim(-0.01, 0.03)

        ax2.set_xlabel(lbls["x"])
        ax1.set_title(r"$\mathrm{Potential}\ \Phi$", loc="left")
        ax2.set_title(r"$\mathrm{Sum:\ }(\Phi+\Psi)$", loc="left")

        ax2.minorticks_on()
        lines=[line1, line2, line3]
        self.set_k_lables(fig, lines)

        y1min, y1max = ax1.get_ylim()
        y2min, y2max = ax2.get_ylim()

        if self.x_rec is not None:
            ax1.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
            ax1.text(self.x_rec-0.1, y1max - 0.06*(y1max-y1min), r"$x_\mathrm{rec}$", fontsize=30)
            ax2.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
            ax2.text(self.x_rec-0.1, y2max - 0.06*(y2max-y2min), r"$x_\mathrm{rec}$", fontsize=30)
        
        if self.x_RM is not None:
            ax1.axvline(self.x_RM, ymax=0.925, color="black", ls="dashdot", lw=self.MarkLines)
            ax1.text(self.x_RM-0.1, y1max - 0.06*(y1max-y1min), r"$x_\mathrm{RM}$", fontsize=30)
            ax2.axvline(self.x_RM, ymax=0.925, color="black", ls="dashdot", lw=self.MarkLines)
            ax2.text(self.x_RM-0.1, y2max - 0.06*(y2max-y2min), r"$x_\mathrm{RM}$", fontsize=30)

        save_push(fig, "potentials")

    def monopole_plot(self)->None:
        fig, ax = plt.subplots()
        line1, = ax.plot(self.k1_x, 4*self.k1_T0, lw=self.LWpoles, color=Colors["k1"], label=lbls["k1"])
        line2, = ax.plot(self.k2_x, 4*self.k2_T0, lw=self.LWpoles, color=Colors["k2"], label=lbls["k2"])
        line3, = ax.plot(self.k3_x, 4*self.k3_T0, lw=self.LWpoles, color=Colors["k3"], label=lbls["k3"])
        # ax.grid(True)
        ax.set_xlim(-15,0)
        ax.set_xlabel(lbls["x"])
        ax.set_title(r"$\mathrm{Overdensity}\ \delta_\gamma=4\Theta_0$", loc="left")
        ax.minorticks_on()
        lines=[line1, line2, line3]
        self.set_k_lables(fig, lines)

        ymin, ymax = ax.get_ylim()
        if self.x_rec is not None:
            ax.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
            ax.text(self.x_rec-0.1, ymax - 0.06*(ymax-ymin), r"$x_\mathrm{rec}$", fontsize=30)
        if self.x_RM is not None:
            ax.axvline(self.x_RM, ymax=0.925, color="black", ls="dashdot", lw=self.MarkLines)
            ax.text(self.x_RM-0.1, ymax - 0.06*(ymax-ymin), r"$x_\mathrm{RM}$", fontsize=30)
        
        save_push(fig, "monopole")

    def dipole_plot(self)->None:
        fig, ax = plt.subplots()
        line1, = ax.plot(self.k1_x, -3*self.k1_T1, lw=self.LWpoles, color=Colors["k1"], label=lbls["k1"])
        line2, = ax.plot(self.k2_x, -3*self.k2_T1, lw=self.LWpoles, color=Colors["k2"], label=lbls["k2"])
        line3, = ax.plot(self.k3_x, -3*self.k3_T1, lw=self.LWpoles, color=Colors["k3"], label=lbls["k3"])
        
        ax.set_xlabel(lbls["x"])
        ax.set_title(r"$\mathrm{Velocity}\ v_\gamma=-3\Theta_1$", loc="left")
        ax.set_xlim(-15,0)
        # ax.grid(True)
        
        ax.minorticks_on()
        lines=[line1, line2, line3]
        self.set_k_lables(fig, lines)
        
        ymin, ymax = ax.get_ylim()
        if self.x_rec is not None:
            ax.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
            ax.text(self.x_rec-0.1, ymax - 0.06*(ymax-ymin), r"$x_\mathrm{rec}$", fontsize=30)
        if self.x_RM is not None:
            ax.axvline(self.x_RM, ymax=0.925, color="black", ls="dashdot", lw=self.MarkLines)
            ax.text(self.x_RM-0.1, ymax - 0.06*(ymax-ymin), r"$x_\mathrm{RM}$", fontsize=30)
        
        save_push(fig, "dipole")

    def quadrapole_plot(self)->None:
        fig, ax = plt.subplots()
        line1, = ax.plot(self.k1_x, self.k1_T2, lw=self.LWpoles, color=Colors["k1"], label=lbls["k1"])
        line2, = ax.plot(self.k2_x, self.k2_T2, lw=self.LWpoles, color=Colors["k2"], label=lbls["k2"])
        line3, = ax.plot(self.k3_x, self.k3_T2, lw=self.LWpoles, color=Colors["k3"], label=lbls["k3"])

        ax.set_xlabel(lbls["x"])
        ax.set_title(r"$\mathrm{Quadrupole}\ \Theta_2$", loc="left")
        ax.set_xlim(-15,0)
        # ax.grid(True)

        ax.minorticks_on()
        lines=[line1, line2, line3]
        self.set_k_lables(fig, lines)
        
        ymin, ymax = ax.get_ylim()
        if self.x_rec is not None:
            ax.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
            ax.text(self.x_rec-0.1, ymax - 0.06*(ymax-ymin), r"$x_\mathrm{rec}$", fontsize=30)
        if self.x_RM is not None:
            ax.axvline(self.x_RM, ymax=0.925, color="black", ls="dashdot", lw=self.MarkLines)
            ax.text(self.x_RM-0.1, ymax - 0.06*(ymax-ymin), r"$x_\mathrm{RM}$", fontsize=30)
        
        save_push(fig, "quadrapole")

    def delta_plot(self)->None:
        fig, ax = plt.subplots(figsize=(12,10))
        line1, = ax.plot(self.k1_x, abs(self.k1_delta_cdm), color=Colors["k1"], label=lbls["k1"])
        line2, = ax.plot(self.k2_x, abs(self.k2_delta_cdm), color=Colors["k2"], label=lbls["k2"])
        line3, = ax.plot(self.k3_x, abs(self.k3_delta_cdm), color=Colors["k3"], label=lbls["k3"])
        ax.plot(self.k1_x, abs(self.k1_delta_b), lw=2, color=Colors["k1"], ls="dashed")
        ax.plot(self.k2_x, abs(self.k2_delta_b), lw=2, color=Colors["k2"], ls="dashed")
        ax.plot(self.k3_x, abs(self.k3_delta_b), lw=2, color=Colors["k3"], ls="dashed")
        ax.set_title(r"$\mathrm{Overdensities\ } \delta_c\ \mathrm{and}\ \delta_b$", loc="left")

        line_lab, = ax.plot(np.nan, np.nan, color="grey")
        dash_lab, = ax.plot(np.nan, np.nan, color="grey", ls="dashed")
        ax.set_xlabel(lbls["x"])
        ax.set_ylim(1e-2, 1e5)
        ax.set_xlim(-15,0)
        ax.set_yscale("log")
        ax.minorticks_on()
        
        lines=[line1, line2, line3]
        self.set_k_lables(fig, lines)
        ax.legend([line_lab, dash_lab], [r"$|\delta_c|$", r"$|\delta_b|$"], fancybox=True, loc="best")

        if self.x_rec is not None:
            ax.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
        if self.x_RM is not None:
            ax.axvline(self.x_RM, color="black", ls="dashdot", lw=self.MarkLines)
    


        
        # ax.legend(loc="best", fancybox=True)
        
        save_push(fig, "delta")

    def delta_comp_plot(self)->None:
        fig, ax = plt.subplots(figsize=(12,10))
        ax.plot(self.k3_x, abs(self.k3_delta_cdm), lw=4, color=Colors["k3c1"], label=r"$|\delta_c|$")
        ax.plot(self.k3_x, abs(self.k3_delta_b), lw=3, color=Colors["k3c2"], label=r"$|\delta_b|$", ls="dashed")
        ax.plot(self.k3_x, abs(4*self.k3_T0), lw=2, color=Colors["k3c3"], label=r"$|\delta_\gamma|$", ls="dashdot")

        ax.legend(fancybox=True, loc="best")

        ax.set_title(r"$\mathrm{Overdensities\ when\ } k=0.1 / \mathrm{Mpc}$", loc="left")

        ax.set_xlabel(lbls["x"])
        ax.set_xlim(-15,0)
        ax.set_ylim(1e-2, 1e5)
        ax.set_yscale("log")
        ax.minorticks_on()

        if self.x_rec is not None:
            ax.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
        if self.x_RM is not None:
            ax.axvline(self.x_RM, color="black", ls="dashdot", lw=self.MarkLines)

        


        save_push(fig, "delta_comparison")

    def velocity_plot(self)->None:
        fig, ax = plt.subplots(figsize=(12,10))
        line1, = ax.plot(self.k1_x, abs(self.k1_v_cdm), color=Colors["k1"], label=lbls["k1"])
        line2, = ax.plot(self.k2_x, abs(self.k2_v_cdm), color=Colors["k2"], label=lbls["k2"])
        line3, = ax.plot(self.k3_x, abs(self.k3_v_cdm), color=Colors["k3"], label=lbls["k3"])
        ax.plot(self.k1_x, abs(self.k1_v_b), lw=2, color=Colors["k1"], ls="dashed")
        ax.plot(self.k2_x, abs(self.k2_v_b), lw=2, color=Colors["k2"], ls="dashed")
        ax.plot(self.k3_x, abs(self.k3_v_b), lw=2, color=Colors["k3"], ls="dashed")

        ax.set_title(r"$\mathrm{Velocities\ } v_c\ \mathrm{and}\ v_b$", loc="left")

        ax.set_xlabel(lbls["x"])
        ax.set_xlim(-15,0)
        ax.set_ylim(1e-4, 5e1)
        ax.set_yscale("log")
        ax.minorticks_on()
        # ax.legend(loc="best", fancybox=True)
        line_lab, = ax.plot(np.nan, np.nan, color="grey")
        dash_lab, = ax.plot(np.nan, np.nan, color="grey", ls="dashed")
        lines=[line1, line2, line3]
        self.set_k_lables(fig, lines)
        ax.legend([line_lab, dash_lab], [r"$|v_c|$", r"$|v_b|$"], fancybox=True, loc="best")
        if self.x_rec is not None:
            ax.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
        if self.x_RM is not None:
            ax.axvline(self.x_RM, color="black", ls="dashdot", lw=self.MarkLines)
        
        save_push(fig, "velocity")


    def velocity_comp_plot(self)->None:
        fig, ax = plt.subplots(figsize=(12,10))
        ax.plot(self.k3_x, abs(self.k3_v_cdm), lw=4, color=Colors["k3c1"], label=r"$|v_c|$")
        ax.plot(self.k3_x, abs(self.k3_v_b), lw=3, color=Colors["k3c2"], label=r"$|v_b|$", ls="dashed")
        ax.plot(self.k3_x, abs(-3*self.k3_T1), lw=2, color=Colors["k3c3"], label=r"$|v_\gamma|$", ls="dashdot")

        ax.legend(fancybox=True, loc="best")

        ax.set_title(r"$\mathrm{Velocities\ when\ } k=0.1 / \mathrm{Mpc}$", loc="left")

        ax.set_xlabel(lbls["x"])
        ax.set_xlim(-15,0)
        ax.set_ylim(1e-3, 5e1)
        ax.set_yscale("log")
        ax.minorticks_on()
        if self.x_rec is not None:
            ax.axvspan(self.x_rec-RECTOL, self.x_rec+RECTOL, color="grey", alpha=0.3)
        if self.x_RM is not None:
            ax.axvline(self.x_RM, color="black", ls="dashdot", lw=self.MarkLines)
   


        save_push(fig, "velocity_comparison")


    def plot_integrand_test(self)->None:
        fig, ax = plt.subplots()
        line1, = ax.plot(self.k1_x, self.k1_Sj5, color=Colors["k1"], label=lbls["k1"])
        line2, = ax.plot(self.k2_x, self.k2_Sj5, color=Colors["k2"], label=lbls["k2"])
        line3, = ax.plot(self.k3_x, self.k3_Sj5, color=Colors["k3"], label=lbls["k3"])
        ax.set_xlabel(lbls["x"])
        ax.set_xlim(-15,0)
        ax.minorticks_on()


        lines=[line1, line2, line3]
        # self.set_k_lables(fig, lines)
        

        save_push(fig, "integrand")



    def make_plots(self)->None:
        self.potentials_plot()
        self.monopole_plot()
        self.dipole_plot()
        self.quadrapole_plot()
        self.delta_plot()
        self.delta_comp_plot()
        self.velocity_plot()
        self.velocity_comp_plot()

        # self.plot_integrand_test()

if __name__=="__main__":
    Pert = Perturbation(["perturbations_k0.001.csv", "perturbations_k0.01.csv", "perturbations_k0.1.csv"])
    Pert.make_plots()
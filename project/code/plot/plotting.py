from plot_utils import *
import analysisM1 as M1
import analysisM2 as M2
import analysisM3 as M3

def runM1():
    M1.testing_Omegas()
    M1.testing_Hp()
    M1.testing_eta()
    M1.conformal_hubble_factor()
    M1.cosmic_conformal_time()
    M1.supernova_data()
    M1.prob_plots()
    M1.goodness_of_fit()
    M1.omega_restrictions_plot()
    M1.posterior_pdf()
    # M1.create_table()

def runM2():
    Rec = M2.Recombination(Data("recombination.csv"))
    Rec.make_plots()
    # Rec.make_table()
    # print(Rec.AnalData)

def runM3():
    Rec = M2.Recombination(Data("recombination.csv"))
    Pert = M3.Perturbation(["perturbations_k0.001.csv", "perturbations_k0.01.csv", "perturbations_k0.1.csv"], Rec.x_rec, M1.x_RM)
    Pert.make_plots()

if __name__ == "__main__":
    runM3()

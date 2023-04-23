#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  h = 0.7;
  OmegaB = 0.05;
  OmegaCDM = 0.45;
  Neff = 0;

  // Recombination parameters
  // double Yp          = 0.245;
  double Yp          = 0;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  // BackgroundCosmology cosmo(.7, .25, .25, 0, Neff, TCMB);

  Utils::StartTiming("Solve & write");
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  // cosmo.output("data/backgroundcosmology.csv");
  // cosmo.output("data/backgroundcosmologyLumDist.csv", log(0.01+1)-1, log(1.3+1)+1, (int)1e5);
  // cosmo.write_table_of_important_values("data/table_of_values.csv");
  Utils::EndTiming("Solve & write");

  // BackgroundCosmology bestFit(0.702, 0.05, 0.259-0.05, 0.067, Neff, TCMB);
  // Utils::StartTiming("Solve best params");
  // bestFit.solve();
  // bestFit.info();
  // bestFit.output("data/bestFitBackground.csv", log(0.01+1)-1, log(1.3+1)+1, (int)1e5);
  // Utils::EndTiming("Solve best params");


  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h

  // Utils::StartTiming("SupernovaFit");
  // mcmc_fit_to_supernova_data("data/supernovadata.txt", "data/results_supernovafitting.csv", "data/best_fit_params.csv");
  // Utils::EndTiming("SupernovaFit");

  // Remove when module is completed
  // return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  // rec.output("data/recombination.csv");
  // rec.analysis_output("data/recomb_analysis.csv");
  
  // Remove when module is completed
  // return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "data/perturbations_k0.01.csv");
  kvalue = 0.1 / Constants.Mpc;
  pert.output(kvalue, "data/perturbations_k0.1.csv");
  kvalue = 0.001 / Constants.Mpc;
  pert.output(kvalue, "data/perturbations_k0.001.csv");
  
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}

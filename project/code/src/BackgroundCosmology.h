#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "Utils.h"

using Vector = std::vector<double>;

class BackgroundCosmology{
  private:
   
    // Cosmological parameters
    double h;                       // Little h = H0/(100km/s/Mpc)
    double OmegaB;                  // Baryon density today
    double OmegaCDM;                // CDM density today
    double OmegaLambda;             // Dark energy density today
    double Neff;                    // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double TCMB;                    // Temperature of the CMB today in Kelvin
   
    // Derived parameters
    double OmegaR;                  // Photon density today (follows from TCMB)
    double OmegaNu;                 // Neutrino density today (follows from TCMB and Neff)
    double OmegaK;                  // Curvature density = 1 - OmegaM - OmegaR - OmegaNu - OmegaLambda
    double OmegaM;                  // Total matter density
    double OmegaR_tot;                  // Total radiation density

    

    // Present values
    double t0;                      // Age of universe today [s] 
    double eta0;                    // Conformal time today [s]

    // Misc
    double toGyr = 1/(1e9*365*24*60*60);  // Conversion from [s] to [Gyr]

    // Start and end of x-integration (can be changed)
    double x_start = Constants.x_start;
    double x_end   = Constants.x_end;

    // Splines to be made
    Spline eta_of_x_spline{"eta"};
    Spline t_of_x_spline{"t"};
 
    //  Private function for easier implementation
    double u_of_x(double x) const;
    double dudx_of_x(double x) const;
    double dduddx_of_x(double x) const;

    // alpha_values
    // struct{
    //   double M=1.0;
    //   double k=0.0;
    //   double R=2.0;
    //   double L=-2.0;
    // }alpha;
    
  public:

    // Constructors 
    BackgroundCosmology() = delete;
    BackgroundCosmology(
        double h, 
        double OmegaB, 
        double OmegaCDM, 
        double OmegaK,
        double Neff, 
        double TCMB
        );

    double H0;                      // The Hubble parameter today H0 = 100h km/s/Mpc
    // Print some useful info about the class

    // Equalities
    double x_RM;                     // Radiation-matter equality
    double x_ML;                     // Matter-dark energy equality
    double x_acc;                    // Acceleration starts
    double x0 = 0.0;                 // Value today
    
    void info() const;

    // Do all the solving
    void solve(int nr_points = (int)1e5);

    // Output some results to file
    void output(const std::string filename, double x_min=-20.0, double x_max = 5.0, int n_pts = (int)1e5) const;

    void write_table_of_important_values(const std::string filename) const;

    // Get functions
    double eta_of_x(double x) const;
    double deta_of_x(double x) const;
    double H_of_x(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;
    // double Hp_ratio(double x) const;
    double get_OmegaB(double x = 0.0) const; 
    double get_OmegaM(double x = 0.0) const; 
    double get_OmegaR(double x = 0.0) const;
    double get_OmegaRtot(double x = 0.0) const; 
    double get_OmegaNu(double x = 0.0) const;
    double get_OmegaCDM(double x = 0.0) const; 
    double get_OmegaLambda(double x = 0.0) const; 
    double get_OmegaK(double x = 0.0) const; 
    double get_OmegaMnu(double x = 0.0) const; 
    double get_H0() const;
    double get_h() const;
    double get_Neff() const;
    double get_TCMB(double x = 0.0) const;

    // Distance measures
    double get_luminosity_distance_of_x(double x) const;
    double get_comoving_distance_of_x(double x) const;
    double get_angular_distance_of_x(double x) const;

    // Extra
    double get_r_of_x(double x) const;
    double get_t_of_x(double x) const;
    double get_z_of_x(double x) const;
  

};

#endif

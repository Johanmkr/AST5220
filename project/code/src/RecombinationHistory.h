#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "BackgroundCosmology.h"

using Vector = std::vector<double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;
 
    // The start and end points for recombination arrays (can be modified)
    const double x_start  = Constants.x_start;
    const double x_end    = Constants.x_end;
    
    // Numbers of points of Xe,ne array (modify as you see fit)
    const int npts_rec_arrays = int(1e6);
  
    // Xe for when to switch between Saha and Peebles
    const double Xe_saha_limit = 0.99;

    // Constants from cosmology class
    double OmegaB; 
    double TCMB;
    double H0;
    double rho_c0;   // Critical density today
    void set_cosmo_constant();

    // // Derived Constants
    // double const_nb_inv;; // Constant part of 1/nb
    // double meTbpow;
    // double const_saha_eq; // Constant part of the saha equation
    // double const_eps_tcmb;  // Constant part of exponent of Saha equation.
    // void set_derived_constants();

    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe 
    void solve_number_density_electrons();
    
    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================

    // The two things we need to solve: Xe/ne and tau
    void solve_for_optical_depth_tau();

    // Splines contained in this class
    Spline log_Xe_of_x_spline{"Xe"};
    Spline log_ne_of_x_spline{"ne"};
    Spline tau_of_x_spline{"tau"}; 
    Spline g_tilde_of_x_spline{"g"};  

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp);
      
    // Utility functions
    // double helper_Saha(double x) const; // Helper function in order to solve the Saha equation

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output some data to file
    void output(const std::string filename) const;


    // Get functions that we must implement
    double tau_of_x(double x) const;
    double dtaudx_of_x(double x) const;
    double ddtauddx_of_x(double x) const;
    double g_tilde_of_x(double x) const;
    double dgdx_tilde_of_x(double x) const;
    double ddgddx_tilde_of_x(double x) const;
    double Xe_of_x(double x) const;
    double ne_of_x(double x) const;
    double get_nb_of_x(double x) const;
    double get_nH_of_x(double x) const;
    double get_Yp() const;
};

#endif

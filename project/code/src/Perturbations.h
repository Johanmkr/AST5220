#ifndef _PERTURBATIONS_HEADER
#define _PERTURBATIONS_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <vector>
#include <fstream>
#include <algorithm>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class Perturbations{
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
   
    // The scales we integrate over
    const int n_k        = 100;
    const double k_min   = Constants.k_min;
    const double k_max   = Constants.k_max;
    
    // Start and end of the time-integration
    const int n_x        = 10000;
    const double x_start = Constants.x_start;
    const double x_end   = Constants.x_end;

    // Below is a full list of splines you probably need, 
    // but you only need to make the splines you will need

    // Splines of scalar perturbations quantities
    Spline2D delta_cdm_spline{"delta_cdm_spline"};
    Spline2D delta_b_spline{"delta_b_spline"};
    Spline2D v_cdm_spline{"v_cdm_spline"};
    Spline2D v_b_spline{"v_b_spline"};
    Spline2D Phi_spline{"Phi_spline"};
    Spline2D Pi_spline{"Pi_spline"};
    Spline2D Psi_spline{"Psi_spline"};
   
    // Splines of source functions (ST for temperature; SE for polarization)
    Spline2D ST_spline{"ST_spline"};
    Spline2D SE_spline{"SE_spline"};
    
    // Splines of mulipole quantities
    // NB: If you use there you have to allocate the container first
    // e.g. Theta_spline = std::vector<Spline2D>(n_ell_theta); before using it
    std::vector<Spline2D> Theta_spline;

    
    //==========================================================
    // [1] Tight coupling ODE system
    //==========================================================

    /**
     * @brief Set the initial conditions at the start (which is in tight coupling)
     * 
     * @param x time
     * @param k k-mode
     * @return Vector 
     */
    Vector set_ic(const double x, const double k) const;
    
    /**
     * @brief Compute the right hand side of the ODE in the tight coupling regime.
     * 
     * @param x time
     * @param k k-mode
     * @param y solution array
     * @param dydx rhs to be set
     * @return int 
     */
    int rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx);
    

    /**
     * @brief Compute time when tight-coupling ends
     * 
     * @param k k-mode
     * @param x_arr x-array
     * @return int 
     */
    int get_tight_coupling_time_idx(const double k, Vector x_arr) const;
    
    //==========================================================
    // [2] The full ODE system 
    //==========================================================
    
    /**
     * @brief Set initial condition after tight coupling
     * 
     * @param y_tight_coupling solution array from tight coupling
     * @param x time
     * @param k k-mode
     * @return Vector 
     */
    Vector set_ic_after_tight_coupling(const Vector &y_tight_coupling, const double x, const double k) const;

    /**
     * @brief Right hand side of the ODE in the full regime
     * 
     * @param x time
     * @param k k-mode
     * @param y solution array
     * @param dydx rhs to be set
     * @return int 
     */
    int rhs_full_ode(double x, double k, const double *y, double *dydx);
    
    //==========================================================
    // [3] Integrate the full system
    //==========================================================
    
    /**
     * @brief Integrate perturbations and spline the result
     * 
     */
    void integrate_perturbations();
    
    //==========================================================
    // [4] Compute source functions from the result
    //==========================================================
    
    /**
     * @brief Compute source functions and spline the result
     * 
     */
    void compute_source_functions(bool SW=true, bool ISW=true, bool DOP=true, bool POL=true);

  public:

    // Constructors
    Perturbations() = default;
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec); 

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output info to file
    void output(const double k, const std::string filename) const;

    // Get the quantities we have integrated
    double get_delta_cdm(const double x, const double k) const;
    double get_delta_b(const double x, const double k) const;
    double get_v_cdm(const double x, const double k) const;
    double get_v_b(const double x, const double k) const;
    double get_Phi(const double x, const double k) const;
    double get_Psi(const double x, const double k) const;
    double get_Pi(const double x, const double k) const;
    double get_Theta(const double x, const double k, const int ell) const;

    // Get the derivatives needed
    double get_dv_b(const double x, const double k) const;
    double get_dPhi(const double x, const double k) const;
    double get_dPsi(const double x, const double k) const;
    double get_dTheta(const double x, const double k, const int ell) const;

    // Get the photon source function
    double get_Source_T(const double x, const double k) const;
    // double get_Source_E(const double x, const double k) const;
};

#endif

#include"Perturbations.h"

/*
  CONSTRUCTOR
*/

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

/*
  DO ALL THE SOLVING
*/

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}


void Perturbations::integrate_perturbations(){
    Utils::StartTiming("integrateperturbation");
    /* 
      SET UP k-ARRAYS AND FULL x-ARRAY AND CREATE THE f-ARRAYS TO STORE THE FULL SOLUTIONS UPON SPLINING
    */
    double log_k_min = log(k_min);
    double log_k_max = log(k_max);
    Vector log_k_array = Utils::linspace(log_k_min, log_k_max, n_k); // Linearly spaced logarithmic values
    Vector k_array = exp(log_k_array);

    Vector x_array = Utils::linspace(x_start, x_end, n_x);

    // Declare vectors to store calculated values upon splining:
    Vector f_delta_cdm(n_k*n_x);
    Vector f_delta_b(n_k*n_x);
    Vector f_v_cdm(n_k*n_x);
    Vector f_v_b(n_k*n_x);
    Vector f_Phi(n_k*n_x);
    Vector f_Psi(n_k*n_x);
    std::vector<Vector> f_Theta(Constants.n_ell_theta, Vector(n_x*n_k));

    // Derivative vector:
    // Vector 


    // Loop over all wavenumbers
    #pragma omp parallel for schedule(dynamic, 1) // Dynamic parallelisation across the k-values
    for(int ik = 0; ik < n_k; ik++){

      /*
        FOR EACH k, FIND TIGHT COUPLING TIME, AND CREATE THE NECESSARY x-ARRAYS
      */

      // Progress bar...
      if( (10*ik) / n_k != (10*ik+10) / n_k ) {
        std::cout << (100*ik+100)/n_k << "% " << std::flush;
        if(ik == n_k-1) std::cout << std::endl;
      }

      // Current value of k
      double k = k_array[ik];

      // Find value to integrate to
      int end_tight_idx = get_tight_coupling_time_idx(k, x_array);

      // Useful numbers
      int nr_points_in_tc_regime = end_tight_idx+1;
      int nr_points_in_full_regime = n_x - nr_points_in_tc_regime;
      int first_full_regime_idx = nr_points_in_tc_regime;

      double x_end_tight = x_array [end_tight_idx];

      // Make x_vector for tight_coupling
      Vector x_tc_arr(nr_points_in_tc_regime);
      for(int i=0; i<nr_points_in_tc_regime; i++){
        x_tc_arr[i] = x_array[i];
      }

      // Make x_vector for full system - need one point overlap
      Vector x_full(nr_points_in_full_regime+1);
      for (int i=0; i<=nr_points_in_full_regime; i++){
        x_full[i] = x_array[i+nr_points_in_tc_regime-1];
      }

      /*
        SOLVE TIGHT-COUPLING EQUATIONS
      */

      // Set up initial conditions in the tight coupling regime
      auto y_tight_coupling_ini = set_ic(x_start, k);

      // The tight coupling ODE system
      ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
        return rhs_tight_coupling_ode(x, k, y, dydx);
      };

      ODESolver ode_tc;
      ode_tc.solve(dydx_tight_coupling, x_tc_arr, y_tight_coupling_ini);

      auto solution_tc = ode_tc.get_data();
      auto derivatives_tc = ode_tc.get_derivative_data();
      
      /*
        SOLVE FULL SYSTEM
      */

      // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
      Vector y_tight_coupling = solution_tc[end_tight_idx];
      auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

      // The full ODE system
      ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
        return rhs_full_ode(x, k, y, dydx);
      };

      ODESolver ode_full;
      ode_full.solve(dydx_full, x_full, y_full_ini);

      auto solution_full = ode_full.get_data();
      auto derivative_full = ode_full.get_derivative_data();
      

      /*
        SET INDICES AND FILL f-ARRAYS WITH VALUES FROM BOTH TIGHT-COUPLING AND FULL REGIME
      */

      int delta_cdm_idx = Constants.ind_deltacdm;
      int delta_b_idx = Constants.ind_deltab;
      int v_cdm_idx = Constants.ind_vcdm;
      int v_b_idx = Constants.ind_vb;
      int Phi_idx = Constants.ind_Phi;
      int start_Theta_idx = Constants.ind_start_theta;

      int n_ell_theta = Constants.n_ell_theta;

      // Set the f-vectors with the corresponding value for each x

      // Fill from tight coupling
      for(int ix=0; ix<nr_points_in_tc_regime; ix++){
        int idx = ix + n_x*ik;
        auto y = solution_tc[ix];
        // References to the quantities we are going to set
        double &delta_cdm       =  y[Constants.ind_deltacdm];
        double &delta_b         =  y[Constants.ind_deltab];
        double &v_cdm           =  y[Constants.ind_vcdm];
        double &v_b             =  y[Constants.ind_vb];
        double &Phi             =  y[Constants.ind_Phi];
        double *Theta           = &y[Constants.ind_start_theta];

        // Calculate some relevant quantities
        double x = x_tc_arr[ix];
        double Hp = cosmo->Hp_of_x(x);
        double c = Constants.c;
        double OmegaR = cosmo->get_OmegaR(x);
        double ckHp = c*k/Hp;
        double dtau = rec->dtaudx_of_x(x);
        f_delta_cdm[idx] = delta_cdm;
        f_delta_b[idx] = delta_b;
        f_v_cdm[idx] = v_cdm;
        f_v_b[idx] = v_b;
        f_Phi[idx] = Phi;
        f_Psi[idx] = -Phi - 12.*Hp*Hp / (c*c*k*k) * OmegaR * Theta[2];
        f_Theta[0][idx] = Theta[0];
        f_Theta[1][idx] = Theta[1];
        f_Theta[2][idx] = -20./45.*ckHp / dtau * Theta[1];
        for(int l=3; l<n_ell_theta; l++){
          f_Theta[l][idx] = -l/(2*l+1)*ckHp/dtau * Theta[l-1];
        }
      }

      // Fill from full regime
      for(int ix=nr_points_in_tc_regime; ix<n_x; ix++){
        int idx = ix + n_x*ik;
        int act_ix = ix - nr_points_in_tc_regime + 1; // +1 due to overlap
        auto y = solution_full[act_ix];
        // References to the quantities we are going to set
        double &delta_cdm       =  y[Constants.ind_deltacdm];
        double &delta_b         =  y[Constants.ind_deltab];
        double &v_cdm           =  y[Constants.ind_vcdm];
        double &v_b             =  y[Constants.ind_vb];
        double &Phi             =  y[Constants.ind_Phi];
        double *Theta           = &y[Constants.ind_start_theta];

        // Calculate some relevant quantities
        double x = x_full[act_ix];
        double Hp = cosmo->Hp_of_x(x);
        double c = Constants.c;
        double OmegaR = cosmo->get_OmegaR(x);
        f_delta_cdm[idx] = delta_cdm;
        f_delta_b[idx] = delta_b;
        f_v_cdm[idx] = v_cdm;
        f_v_b[idx] = v_b;
        f_Phi[idx] = Phi;
        f_Psi[idx] = -Phi - 12.*Hp*Hp / (c*c*k*k) * OmegaR * Theta[2];
        for(int l=0; l<n_ell_theta; l++){
          f_Theta[l][idx] = Theta[l];
        }
      }
    }

    Utils::EndTiming("integrateperturbation");

    /*
      MAKE ALL NECESSARY SPLINES
    */
    
    delta_cdm_spline.create(x_array, k_array, f_delta_cdm);
    delta_b_spline.create(x_array, k_array, f_delta_b);
    v_cdm_spline.create(x_array, k_array, f_v_cdm);
    v_b_spline.create(x_array, k_array, f_v_b);
    Phi_spline.create(x_array, k_array, f_Phi);
    Psi_spline.create(x_array, k_array, f_Psi);
    Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
    for(int l=0; l<Constants.n_ell_theta; l++){
      Theta_spline[l].create(x_array, k_array, f_Theta[l]);
    }
  }

/*
  SET INITIAL CONDITION AT THE START (TIGHT COUPLING)
*/

Vector Perturbations::set_ic(const double x, const double k) const{
  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];

  // Set the initial conditions
  double Psi      = -2./3.;
  double c        = Constants.c;
  double Hp       = cosmo->Hp_of_x(x);
  double ckHp     = c*k/Hp;

  // Compute the scalar quantities
  Phi             = -Psi;
  delta_cdm       = -3./2.*Psi;
  delta_b         = -3./2.*Psi;
  v_cdm           = -.5*ckHp*Psi;
  v_b             = -.5*ckHp*Psi;
  // Compute the photon temperature perturbation
  Theta[0]        = -.5*Psi;
  Theta[1]        = ckHp/6.*Psi;

  return y_tc;
}

/*
  SET INITIAL CONDITIONS FOR FULL SYSTEM (AFTER TIGHT COUPLING)
*/

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];

  // Useful quantities
  double c = Constants.c;
  double Hp = cosmo->Hp_of_x(x);
  double ckHp = c*k/Hp;

  // Set the scalar quantities
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  // Set the photon temperature multipoles
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = -20./45.*ckHp/rec->dtaudx_of_x(x)*Theta[1];
  for(int l=3; l<n_ell_theta; l++){
    Theta[l] = -l/(2.*l+1)*ckHp/rec->dtaudx_of_x(x)*Theta[l-1];
  }
  return y;
}

/*
  COMPUTE THE TIME FOR WHICH TIGHT COUPLING ENDS
*/

int Perturbations::get_tight_coupling_time_idx(const double k ,Vector x_arr) const{
  // Initial values:
  int idx_tc = 0;
  bool tight_coupling = true;
  double x_current;
  double ckHp;
  double dtaudx;
  while(tight_coupling){
    x_current = x_arr[idx_tc];
    // Calculate necessary stuff
    ckHp = Constants.c*k / cosmo->Hp_of_x(x_current);
    dtaudx = abs(rec->dtaudx_of_x(x_current));
    // Check each condition
    if(dtaudx<10){
      tight_coupling = false;
    }
    else if(dtaudx<10*ckHp){
      tight_coupling = false;
    }
    else if(x_current>-8.3){
      tight_coupling = false;
    }
    else{
      idx_tc++;
    }
  }

  return idx_tc;
}


/*
  COMPUTE THE PHOTON TEMPERATURE SOURCE FUNCTION
*/

void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  /*
    SET UP x AND k-ARRAYS
  */
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, n_k); // Linearly spaced logarithmic values
  Vector k_array = exp(log_k_array); 

  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector ST_SW_array(k_array.size() * x_array.size());
  Vector ST_ISW_array(k_array.size() * x_array.size());
  Vector ST_DOP_array(k_array.size() * x_array.size());
  Vector ST_POL_array(k_array.size() * x_array.size());


  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int idx             = ix + n_x * ik;

      // Background parameters
      const double Hp           = cosmo->Hp_of_x(x);
      const double dHp          = cosmo->dHpdx_of_x(x);
      const double ddHp         = cosmo->ddHpddx_of_x(x);

      //  Recombination parameters
      const double tau          = rec->tau_of_x(x);
      const double dtau         = rec->dtaudx_of_x(x);
      const double ddtau        = rec->ddtauddx_of_x(x);
      const double g            = rec->g_tilde_of_x(x);
      const double dg           = rec->dgdx_tilde_of_x(x);
      const double ddg          = rec->ddgddx_tilde_of_x(x);

      //  Derivative parameters
      double dPsi               = get_dPsi(x,k);
      double dPhi               = get_dPhi(x,k);
      double dv_b               = get_dv_b(x,k);
      double dT1                = get_dTheta(x,k,1);
      double dT2                = get_dTheta(x,k,2);
      double dT3                = get_dTheta(x,k,3);

      // Other parameters
      double T0                 = get_Theta(x,k,0);
      double T1                 = get_Theta(x,k,1);
      double T2                 = get_Theta(x,k,2);
      double T3                 = get_Theta(x,k,3);
      double v_b                = get_v_b(x,k);
      double Psi                = get_Psi(x,k);


      // Constants
      double c                  = Constants.c;
      double ckHp               = c*k/Hp;

      // Compute the second derivative of the quadrupole
      double ddT2 = 2./5.* ckHp * (dT1 - dHp/Hp*T1)
                  + 3./5.*ckHp * (dHp/Hp*T3-dT3)
                  + 9./10. * (ddtau*T2 + dtau*dT2);
      // double ddT2 = get_ddTheta(x, k, 2);

      //  Calculate terms, starting from 0
      double source_func_val = 0.0;

      // //  + contribution from Sachs-Wolfe (SW) effect
      // if(SW){
      //   source_func_val += g*(T0 + Psi + T2/4.);
      // }

      // //  + contributino from Integrated Sachs-Wolfe (ISW) effect
      // if(ISW){
      //   source_func_val += exp(-tau)*(dPsi-dPhi);
      // }

      // //  + contribution from doppler effect
      // if(DOP){
      //   double DOP_term  = -1./(c*k) *
      //                   (
      //                     Hp*g*dv_b + Hp*dg*v_b + dHp*g*v_b
      //                   );
      //   source_func_val += DOP_term;
      // }

      // //  + contribution from polarisation
      // if(POL){
      //   double POL_term  = 3./(4.*c*c*k*k) *
      //                   (
      //                     g*T2*(ddHp*Hp + dHp*dHp)
      //                     + Hp*Hp*(g*ddT2 + ddg*T2 + 2*dg*dT2)
      //                     + 3*dHp*Hp*(g*dT2 + dg*T2)
      //                   );
      //   source_func_val += POL_term;
      // }

      // ST_array[idx] = source_func_val;




      double SW_term      = g*(T0 + Psi + T2/4.);
      double ISW_term     = exp(-tau)*(dPsi-dPhi);
      double DOP_term     = -1./(c*k) *
                        (
                          Hp*g*dv_b + Hp*dg*v_b + dHp*g*v_b
                        );
      double POL_term     = 3./(4.*c*c*k*k) *
                        (
                          g*T2*(ddHp*Hp + dHp*dHp)
                          + Hp*Hp*(g*ddT2 + ddg*T2 + 2*dg*dT2)
                          + 3*dHp*Hp*(g*dT2 + dg*T2)
                        );
      // Fill corresponding arrays

      ST_SW_array[idx] = SW_term;
      ST_ISW_array[idx] = ISW_term;
      ST_DOP_array[idx] = DOP_term;
      ST_POL_array[idx] = POL_term;
      ST_array[idx] = SW_term + ISW_term + DOP_term + POL_term;

    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  ST_SW_spline.create(x_array, k_array, ST_SW_array, "SW source_temp_x_k");
  ST_ISW_spline.create(x_array, k_array, ST_ISW_array, "ISW source_temp_x_k");
  ST_DOP_spline.create(x_array, k_array, ST_DOP_array, "DOP source_temp_x_k");
  ST_POL_spline.create(x_array, k_array, ST_POL_array, "POL source_temp_x_k");

  Utils::EndTiming("source");
}

/*
  COMPUTE THE RIGHT HAND SIDE OF THE ODE IN THE TIGHT COUPLING REGIME
*/

int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  // const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  // const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  // double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //  Background parameters
  double Hp         = cosmo->Hp_of_x(x);
  double dHp        = cosmo->dHpdx_of_x(x);
  double H0         = cosmo->H0;
  double OmegaCDM   = cosmo->get_OmegaCDM(x);
  double OmegaB     = cosmo->get_OmegaB(x);
  double OmegaR     = cosmo->get_OmegaR(x);

  //  Recombination parameters
  double Rinv       = 1. / rec->get_R_of_x(x);
  double dtau       = rec->dtaudx_of_x(x);
  double ddtau      = rec->ddtauddx_of_x(x);

  //  Useful quantities
  double ckHp       = Constants.c*k/Hp;
  double Theta2     = -20./(45.*dtau)*ckHp*Theta[1];
  double Y          = OmegaCDM*delta_cdm + OmegaB*delta_b + 4.*OmegaR * Theta[0];
  double Psi        = -Phi - 12.*Hp*Hp / (Constants.c*Constants.c * k *k) * OmegaR*Theta2; 

  
  dPhidx            = Psi - 1./3. * ckHp*ckHp*Phi + Y/2.;
  
  dThetadx[0]       = - ckHp * Theta[1] - dPhidx;
  
  ddelta_cdmdx      = ckHp * v_cdm - 3.*dPhidx;
  
  ddelta_bdx        = ckHp * v_b - 3.*dPhidx;
  
  dv_cdmdx          = - v_cdm - ckHp*Psi;

  // Tight coupling equations
  double q          = (
    -(ddtau * (1.+Rinv) + (1.-Rinv)*dtau)*(3.*Theta[1]+v_b)
    - ckHp*Psi + (1.-dHp/Hp)*ckHp*(-Theta[0]+2.*Theta2)-ckHp*dThetadx[0]
  )/
  ((1+Rinv)*dtau + dHp/Hp-1);

  dv_bdx            = (
    - v_b - ckHp*Psi + Rinv * (q + ckHp*(-Theta[0] + 2.*Theta2) - ckHp*Psi)
  )/(1.+Rinv);

  dThetadx[1]       = 1./3.*(q - dv_bdx);

  return GSL_SUCCESS;
}

/*
  COMPUTE THE RIGHT HAND SIDE OF THE ODE IN THE FULL SYSTEM REGIME
*/

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];

  
  //  Background parameters
  double Hp         = cosmo->Hp_of_x(x);
  double dHp        = cosmo->dHpdx_of_x(x);
  double H0         = cosmo->H0;
  double OmegaCDM   = cosmo->get_OmegaCDM(x);
  double OmegaB     = cosmo->get_OmegaB(x);
  double OmegaR     = cosmo->get_OmegaR(x);

  //  Recombination parameters
  double Rinv       = 1. / rec->get_R_of_x(x);
  double dtau       = rec->dtaudx_of_x(x);
  double ddtau      = rec->ddtauddx_of_x(x);

  //  Useful quantities
  double ckHp       = Constants.c*k/Hp;
  double Y          = OmegaCDM*delta_cdm + OmegaB*delta_b + 4.*OmegaR * Theta[0];
  double Psi        = -Phi - 12.*Hp*Hp / (Constants.c*Constants.c * k *k) * OmegaR*Theta[2]; 

  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx          = Psi - 1./3. * ckHp*ckHp*Phi + Y/2.;
  ddelta_cdmdx    = ckHp * v_cdm - 3.*dPhidx;
  ddelta_bdx      = ckHp * v_b - 3.*dPhidx;
  dv_cdmdx        = -v_cdm - ckHp*Psi;
  dv_bdx          = -v_b - ckHp*Psi + dtau*Rinv * (3*Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0]     = -ckHp * Theta[1] - dPhidx;
  dThetadx[1]     = 1./3.*ckHp*Theta[0] - 2./3.*ckHp*Theta[2] + 1./3.*ckHp*Psi + dtau*(Theta[1]+1./3.*v_b);

  int l = 2;
  dThetadx[l]     = l*ckHp*Theta[l-1]/(2.*l+1.) - (l+1.)*ckHp*Theta[l+1]/(2.*l+1.) + dtau * (Theta[l] - Theta[2]/10.);

  for(l=3; l<n_ell_theta-1; l++){
    dThetadx[l]   = l*ckHp*Theta[l-1]/(2.*l+1.) - (l+1.)*ckHp*Theta[l+1]/(2.*l+1.) + dtau * Theta[l];
  }

  l = n_ell_theta-1;
  dThetadx[l]     = ckHp*Theta[l-1] - Constants.c*(l+1.)*Theta[l]/(Hp*cosmo->eta_of_x(x)) + dtau*Theta[l];

  return GSL_SUCCESS;
}

/*
  GET-METHODS
*/

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Theta_spline[2](x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_T_SW(const double x, const double k) const{
  return ST_SW_spline(x,k);
}
double Perturbations::get_Source_T_ISW(const double x, const double k) const{
  return ST_ISW_spline(x,k);
}
double Perturbations::get_Source_T_DOP(const double x, const double k) const{
  return ST_DOP_spline(x,k);
}
double Perturbations::get_Source_T_POL(const double x, const double k) const{
  return ST_POL_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}

/*
  DERIVATIVE GET-METHODS
*/

double Perturbations::get_dv_b(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}
double Perturbations::get_dPhi(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}
double Perturbations::get_dPsi(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}
double Perturbations::get_dTheta(const double x, const double k, const int ell) const{
  return Theta_spline[ell].deriv_x(x,k);
}
double Perturbations::get_ddTheta(const double x, const double k, const int ell) const{
  return Theta_spline[ell].deriv_xx(x,k);
}


/*
  PRINT USEFUL INFO
*/

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

/*
  OUTPUT TO FILE FOR GIVEN k
*/

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  fp << " x , " << " T0 , " << " T1 , " << " T2 , " << " Phi , " << " Psi , "<< " delta_cdm , " << " delta_b , " << " v_cdm , " << " v_b , " << " Pi , " << " S , " << " Sj5 , " << " Sj50 , " << " Sj500  , "<< "\n";
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " , ";
    fp << get_Theta(x,k,0)   << " , ";
    fp << get_Theta(x,k,1)   << " , ";
    fp << get_Theta(x,k,2)   << " , ";
    fp << get_Phi(x,k)       << " , ";
    fp << get_Psi(x,k)       << " , ";
    fp << get_delta_cdm(x,k) << " , ";
    fp << get_delta_b(x,k)   << " , ";
    fp << get_v_cdm(x,k)     << " , ";
    fp << get_v_b(x,k)       << " , ";
    fp << get_Pi(x,k)        << " , ";
    fp << get_Source_T(x,k)  << " , ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " , ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " , ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " , ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}


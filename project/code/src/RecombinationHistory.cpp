#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{} 

//====================================================
// Do all the solving we need to do
//====================================================

// Declare physical constants
const double k_b         = Constants.k_b;
const double G           = Constants.G;
const double m_e         = Constants.m_e;
const double hbar        = Constants.hbar;
const double m_H         = Constants.m_H;
const double epsilon_0   = Constants.epsilon_0;
const double H0_over_h   = Constants.H0_over_h;

const double c           = Constants.c;
const double sigma_T     = Constants.sigma_T;
const double lambda_2s1s = Constants.lambda_2s1s;





// Misc constants
const double global_tol = 1e-4;

void RecombinationHistory::set_cosmo_constant(){
  OmegaB = cosmo->get_OmegaB();
  TCMB = cosmo->get_TCMB();
  H0 = cosmo->get_H0();
  rho_c0 = 3.0*H0 * H0 / (8.0 * M_PI * G);   // Critical density today
}

// void RecombinationHistory::set_derived_constants(){
//   const_nb_inv = m_H / (OmegaB*rho_c0); // Constant part of 1/nb
//   meTbpow = std::pow(k_b*m_e*TCMB/(2.0*M_PI*hbar*hbar), 3.0/2.0); 
//   const_saha_eq = const_nb_inv * meTbpow; // Constant part of the saha equation
//   const_eps_tcmb = epsilon_0 / (TCMB * k_b); // Constant part of exponent of Saha equation.
// }

void RecombinationHistory::solve(){

  // Set constants
  set_cosmo_constant();
  // set_derived_constants();    
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  // DONE?
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  int idx = 0;
  // for(int i = 0; i < npts_rec_arrays; i++){
  while(saha_regime){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[idx]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Fill Xe and ne arrays with current values
    Xe_arr[idx] = Xe_current;
    ne_arr[idx] = ne_current;

    // Update the index count
    idx += 1;
  
    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit || idx>=npts_rec_arrays)
      saha_regime = false;
  }



  // Find number of elements of origianl x_array filled by the Saha equataion (which are valid)
  const int nr_elements_filled_by_Saha = idx-1;
  int last_Saha_idx = idx-2;

  std::cout << "Idx: " << idx << std::endl;
  std::cout << "idx - 1: " << Xe_arr[idx-1] << std::endl;

  std::cout << "last_saha_idx-1 " << Xe_arr[last_Saha_idx-1] << std::endl;
  std::cout << "last_saha_idx-2 " << Xe_arr[last_Saha_idx-2] << std::endl;
  std::cout << "last_saha_idx-3 " << Xe_arr[last_Saha_idx-3] << std::endl;
  std::cout << "last_saha_idx-4 " << Xe_arr[last_Saha_idx-4] << std::endl;
  std::cout << "last_saha_idx-5 " << Xe_arr[last_Saha_idx-5] << std::endl;
  std::cout << "last_saha_idx-6 " << Xe_arr[last_Saha_idx-6] << std::endl;
  std::cout << "last_saha_idx-7 " << Xe_arr[last_Saha_idx-7] << std::endl;
  std::cout << "last_saha_idx-8 " << Xe_arr[last_Saha_idx-8] << std::endl;
  std::cout << "last_saha_idx-9 " << Xe_arr[last_Saha_idx-9] << std::endl;
  std::cout << "last_saha_idx-10 " << Xe_arr[last_Saha_idx-10] << std::endl;




  //  Check if we have already filled the entire array or if we have to solve Peebles.
  if(nr_elements_filled_by_Saha < npts_rec_arrays){

    // Make new x-array of the x-points not solved by Saha.
    const int npts_ode_array = npts_rec_arrays-last_Saha_idx;
    Vector x_ode(npts_ode_array);

    // Fill the new x-array
    for(int i=0; i<npts_ode_array; i++){
      x_ode[i] = x_array[i+last_Saha_idx];
    }

    // The Peebles ODE equation
    ODESolver peebles_Xe_ode;
    ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
      return rhs_peebles_ode(x, Xe, dXedx);
    };

    // Initial Xe will be the last (valid) added value from the Saha regime. 
    double Xe_init_val = Xe_arr[last_Saha_idx];
    Vector Xe_init_vec{Xe_init_val};

    // Solve the ODE
    peebles_Xe_ode.solve(dXedx, x_ode, Xe_init_vec);

    // Get resultant Xe-array (will be of length (npts_ode_array))
    auto Xe_ode = peebles_Xe_ode.get_data_by_component(0);


    // Update original Xe and ne arrays
    for(int i=last_Saha_idx; i<npts_rec_arrays; i++){
      // Calculate values
      double nH_temp, Xe_temp, ne_temp;
      nH_temp = get_nH_of_x(x_array[i]);
      Xe_temp = Xe_ode[i-last_Saha_idx];
      ne_temp = Xe_temp*nH_temp;

      // Fill arrays
      Xe_arr[i] = Xe_temp;
      ne_arr[i] = ne_temp;
    }

    std::cout << "xe ode 0 " << Xe_ode[0] << std::endl;
    std::cout << "xe ode 1 " << Xe_ode[1] << std::endl;
    std::cout << "xe ode 2 " << Xe_ode[2] << std::endl;
    std::cout << "xe ode 3 " << Xe_ode[3] << std::endl;
    std::cout << "xe ode 4 " << Xe_ode[4] << std::endl;
    std::cout << "xe ode 5 " << Xe_ode[5] << std::endl;
    std::cout << "xe ode 6 " << Xe_ode[6] << std::endl;
  }
  else{
    std::cout << "------ Xe-array completely filled by Saha equation. Are you sure this is correct? ------" << std::endl;
  }



  std::cout << "last_saha_idx+23 " << Xe_arr[last_Saha_idx+23] << std::endl;
  std::cout << "last_saha_idx+22 " << Xe_arr[last_Saha_idx+22] << std::endl;
  std::cout << "last_saha_idx+21 " << Xe_arr[last_Saha_idx+21] << std::endl;
  std::cout << "last_saha_idx+20 " << Xe_arr[last_Saha_idx+20] << std::endl;
  std::cout << "last_saha_idx+19 " << Xe_arr[last_Saha_idx+19] << std::endl;
  std::cout << "last_saha_idx+18 " << Xe_arr[last_Saha_idx+18] << std::endl;
  std::cout << "last_saha_idx+17 " << Xe_arr[last_Saha_idx+17] << std::endl;
  std::cout << "last_saha_idx+16 " << Xe_arr[last_Saha_idx+16] << std::endl;
  std::cout << "last_saha_idx+15 " << Xe_arr[last_Saha_idx+15] << std::endl;
  std::cout << "last_saha_idx+14 " << Xe_arr[last_Saha_idx+14] << std::endl;
  std::cout << "last_saha_idx+13 " << Xe_arr[last_Saha_idx+13] << std::endl;
  std::cout << "last_saha_idx+12 " << Xe_arr[last_Saha_idx+12] << std::endl;
  std::cout << "last_saha_idx+11 " << Xe_arr[last_Saha_idx+11] << std::endl;
  std::cout << "last_saha_idx+10 " << Xe_arr[last_Saha_idx+10] << std::endl;
  std::cout << "last_saha_idx+9 " << Xe_arr[last_Saha_idx+9] << std::endl;
  std::cout << "last_saha_idx+8 " << Xe_arr[last_Saha_idx+8] << std::endl;
  std::cout << "last_saha_idx+7 " << Xe_arr[last_Saha_idx+7] << std::endl;
  std::cout << "last_saha_idx+6 " << Xe_arr[last_Saha_idx+6] << std::endl;
  std::cout << "last_saha_idx+5 " << Xe_arr[last_Saha_idx+5] << std::endl;
  std::cout << "last_saha_idx+4 " << Xe_arr[last_Saha_idx+4] << std::endl;
  std::cout << "last_saha_idx+3 " << Xe_arr[last_Saha_idx+3] << std::endl;
  std::cout << "last_saha_idx+2 " << Xe_arr[last_Saha_idx+2] << std::endl;
  std::cout << "last_saha_idx+1 " << Xe_arr[last_Saha_idx+1] << std::endl;
  std::cout << "last_saha_idx " << Xe_arr[last_Saha_idx] << std::endl;
  std::cout << "last_saha_idx-1 " << Xe_arr[last_Saha_idx-1] << std::endl;
  std::cout << "last_saha_idx-2 " << Xe_arr[last_Saha_idx-2] << std::endl;
  std::cout << "last_saha_idx-3 " << Xe_arr[last_Saha_idx-3] << std::endl;
  std::cout << "last_saha_idx-4 " << Xe_arr[last_Saha_idx-4] << std::endl;
  std::cout << "last_saha_idx-5 " << Xe_arr[last_Saha_idx-5] << std::endl;
  std::cout << "last_saha_idx-6 " << Xe_arr[last_Saha_idx-6] << std::endl;
  std::cout << "last_saha_idx-7 " << Xe_arr[last_Saha_idx-7] << std::endl;
  std::cout << "last_saha_idx-8 " << Xe_arr[last_Saha_idx-8] << std::endl;
  std::cout << "last_saha_idx-9 " << Xe_arr[last_Saha_idx-9] << std::endl;
  std::cout << "last_saha_idx-10 " << Xe_arr[last_Saha_idx-10] << std::endl;



  // Create new arrays containing the log of Xe and ne
  Vector log_Xe_arr(npts_rec_arrays);
  Vector log_ne_arr(npts_rec_arrays);

  for(int i=0; i<npts_rec_arrays; i++){
    log_Xe_arr[i] = log(Xe_arr[i]);
    log_ne_arr[i] = log(ne_arr[i]);
  }

  // Make splines of the result
  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "Xe");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  
  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  // DONE?

  double Tb = cosmo->get_TCMB(x);
  double nb = get_nb_of_x(x);

  double b = 1./nb * std::pow(k_b*m_e*Tb/(2.*M_PI*hbar*hbar), 1.5) * exp(-epsilon_0/Tb);
  std::cout<< b<< std::endl;
  if(b>1e7)
    Xe = 1.0;
  else{
    Xe = (-b+sqrt(b*b+4.0*b))/2.0;
  }

  // double b = const_saha_eq * exp(helper_Saha(x));
  // double discr = b*b+4.0*b;
  // if(discr>1e9){
  //   Xe = 1.0;
  // }
  // else{
  //   Xe = (-b+sqrt(discr))/2;
  // }
  // double root_val;
  // if(abs(discr-1.0)<global_tol){
  //   root_val = 1.0 + (discr-1.0)/2.0;
  // }
  // else{
  //   root_val = sqrt(discr);
  // }
  // Xe = (-b+root_val)/2.0;
  // if(Xe>1e6)
  // std::cout<< b << std::endl;
  // std::cout<< root_val << std::endl;

  const double nH = get_nH_of_x(x);
  ne = Xe*nH;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Constants for finding RHS of peebles eq.
  const double Tb = cosmo->get_TCMB(x);
  const double eps_tb = epsilon_0/Tb;
  const double alpha = c/hbar*sqrt(3.0*sigma_T/(8.0*M_PI))*m_e; // Dimensionless fine structure constant. 
  const double H = cosmo->H_of_x(x);

  // Finding the terms involved in the RHS
  const double phi2 = 0.448 * log(eps_tb); // dimensionless
  const double alpha2 = hbar*hbar/c*64.0*M_PI*alpha*alpha*phi2 / (m_e * m_e) * sqrt(eps_tb/(27.0*M_PI)); // dimension m^3/s
  const double beta = alpha2* std::pow(k_b*m_e*Tb/(2.0*M_PI), 1.5) * exp(-eps_tb);
  double beta2; // dimension 1/s

  // Checking for large exponent to avoid overflow
  if(eps_tb > 200.0){
    beta2 = 0.0;
  }
  else{
    beta2 = beta*exp(eps_tb*3.0/4.0); // dimension 1/s
  }
  // beta2 = beta*exp(eps_tb*3/4); // dimension 1/s

  const double nH = rho_c0*OmegaB / (m_H * a * a * a); // dimension 1/m^3
  const double n1s = (1.0-X_e)*nH; // dimension 1/m^3
  const double Lambda_alpha = 1.0 / (hbar*hbar*hbar*c*c*c) * H * 27.0 * epsilon_0*epsilon_0*epsilon_0 / (64.0 * M_PI * M_PI * n1s); // dimension 1/s
  const double Cr = (lambda_2s1s + Lambda_alpha) / (lambda_2s1s + Lambda_alpha + beta2); // dimensionless

  const double rhs = Cr/H * (beta*(1.0-X_e) - nH*alpha2*X_e*X_e);

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...
  // DONE?
  
  dXedx[0] = rhs;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = int(1e6);

  // BELOW DOES NOT WORK WHY?
  // Vector x_array_1 = Utils::linspace(0,7, npts); //after recombination
  // Vector x_array_2 = Utils::linspace(7,9, 2*npts); //for recombination
  // Vector x_array_3 = Utils::linspace(9,20, npts); // before recombination

  // Vector x_array_tau_rev(4*npts);

  // for(int i=0; i<npts;i++){
  //   x_array_tau_rev[i] = x_array_1[i];
  // }
  // for(int i=npts; i<3*npts; i++){
  //   x_array_tau_rev[i] = x_array_2[i-npts];
  // }
  // for(int i=3*npts; i<4*npts;i++){
  //   x_array_tau_rev[i] = x_array_3[i-3*npts];
  // }

  // x-array that goes backwards in time, by using positive, strictly increasing x-values
  Vector x_array_tau_rev = Utils::linspace(0, -x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODESolver tauOde;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth

    // Derivative for -dtaudx

    // Since we are using positive increasing x-values, we must provide the negative x-values in the functions for ne_of_x and H_of_x. Notice that this is actually -dtaudx in order for ut to integrate backwards. 
    
    dtaudx[0] = c*sigma_T*ne_of_x(-x)/(cosmo->H_of_x(-x));

    return GSL_SUCCESS;
  };

  // Initial Xe will be the last added value from the Saha regime. 
  Vector tau_init_vec{0.0};

  tauOde.solve(dtaudx, x_array_tau_rev, tau_init_vec);

  auto tau_vec_inv = tauOde.get_data_by_component(0);

  // USEFULL ANYMORE?
  // This spline is now for positive x-values -> need to account for this when calling the function
  // tau_of_x_spline.create(x_array_tau_rev, tau_vec_inv, "tau");

  // Reverse the arrays
  Vector x_array_tau(npts);
  Vector tau_vec(npts);
  for(int i=0;i<npts;i++){
    x_array_tau[i] = -x_array_tau_rev[npts-1-i];
    tau_vec[i] = tau_vec_inv[npts-1-i];
  }

  // spline the result
  tau_of_x_spline.create(x_array_tau, tau_vec, "tau");


  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...
  // DONE?

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Vector g_tilde(npts);
  for(int i=0;i<npts;i++){
    double const x_loc = x_array_tau[i];
    g_tilde[i] = -dtaudx_of_x(x_loc) * exp(-tau_of_x(x_loc));
  }

  g_tilde_of_x_spline.create(x_array_tau, g_tilde, "g");



  Utils::EndTiming("opticaldepth");
}

//====================================================
// Utility methods
//====================================================
// double RecombinationHistory::helper_Saha(double x) const{
//   return 3.0/2.0 * x - const_eps_tcmb*exp(x);
// }


//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_nb_of_x(double x) const{
  return OmegaB * rho_c0 / m_H * exp(-3.0*x);
}

double RecombinationHistory::get_nH_of_x(double x) const{
  return get_nb_of_x(x); //no helium
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = int(1e6);
  const double x_min   = -12.0;
  const double x_max   = 0;

  fp << "  x  , " << "   Xe  , " << "  ne   , " << " tau  , " << "  dtaudx    , " << " ddtauddx  , " << "  g    , " << "    dgdx  , " <<  "  ddgddx   ," << "\n";

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " , ";
    fp << Xe_of_x(x)           << " , ";
    fp << ne_of_x(x)           << " , ";
    fp << tau_of_x(x)          << " , ";
    fp << dtaudx_of_x(x)       << " , ";
    fp << ddtauddx_of_x(x)     << " , ";
    fp << g_tilde_of_x(x)      << " , ";
    fp << dgdx_tilde_of_x(x)   << " , ";
    fp << ddgddx_tilde_of_x(x) << " , ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}


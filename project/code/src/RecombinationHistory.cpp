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

// Import constant from cosmology class
const double OmegaB = cosmo->get_OmegaB();
const double TCMB = cosmo->get_TCMB();
const double H0 = cosmo->get_H0();

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


// Derived Constants
const double rho_c0 = 3*H0 * H0 / (8 * M_PI * G);   // Critical density today

// Created constants
const double const_nb_inv = m_H / (OmegaB*rho_c); // Constant part of 1/nb
const double meTbpow = std::pow(m_e*TCMB/(2*M_PI), 3/2);
const double const_saha_eq = const_nb_inv * meTbpow; // Constant part of the saha equation
const double const_eps_tcmb = epsilon_0 / TCMB; // Constant part of exponent of Saha equation.

// Misc constants
const double global_tol = 1e-5;


void RecombinationHistory::solve(){
    
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
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      //...
      //...

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...
    
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...

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
  
  const double b = const_saha_eq * exp(helper_Saha(x));
  const double discr = b*b+4*b;
  const double root_val;
  if(abs(discr-1)<tol){
    root_val = 1 + (discr-1)/2;
  }
  else{
    root_val = sqrt(discr);
  }
  Xe = (-b+root_val)/2;
  const double nH = exp(-3*x)/const_nb_inv;
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
  const double Tb = TCMB / a;
  const double eps_tb = const_eps_tcmb * a;
  const double alpha = sqrt(3*sigma_T/(8*M_PI))*m_e;
  const double H = cosmo->H_of_x(x);

  // Finding the terms involved in the RHS
  const double phi2 = 0.448 * log(eps_tb);
  const double alpha2 = 64*M_PI*alpha*alpha*phi2 / (m_e * m_e) * sqrt(eps_tb/(27*M_PI));
  const double beta = alpha2*meTbpow / (sqrt(a)*sqrt(a)*sqrt(a))*exp(-eps_tb);
  const dobule beta2;
  if(eps_tb > 200){
    beta2 = 0;
  }
  else{
    beta2 = beta*exp(eps_tb*3/4)
  }
  const double nH = 1 / (a*a*a*const_nb_inv);
  const double n1s = (1-X_e)*nH;
  const double Lambda_alpha = H * 9 * epsilon_0*epsilon_0*epsilon_0 / (64 * M_PI * M_PI * n1s);
  const double Cr = (lambda_2s1s + Lambda_alpha) / (lambda_2s1s + Lambda_alpha + beta2);

  const double rhs = Cr/H * (beta*(1-X_e) - nH*alpha2*X_e*X_e);

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
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Utility methods
//====================================================
double RecombinationHistory::helper_Saha(double x) const{
  return 3/2 * x -const_eps_tcmb*exp(x);
}


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

  return 0.0;
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
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

  return 0.0;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
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
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}


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
const double global_tol = 1e-7;

void RecombinationHistory::set_cosmo_constant(){
  OmegaB = cosmo->get_OmegaB();
  TCMB = cosmo->get_TCMB();
  H0 = cosmo->get_H0();
  rho_c0 = 3.0*H0 * H0 / (8.0 * M_PI * G);   // Critical density today
}

void RecombinationHistory::solve(){

  // Set constants
  set_cosmo_constant();
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
  solve_for_sound_horizon();
}

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  // Set up arrays to store Xe and Ne (and XeSaha)
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_arr_saha(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  int idx = 0;
  while(saha_regime){
    //  Get electorn fraction from Saha equation
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[idx]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Fill Xe and ne arrays with current values
    Xe_arr[idx] = Xe_current;
    Xe_arr_saha[idx] = Xe_current;
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

  // Fill rest with saha for comparison
  for(int i=last_Saha_idx+1; i<npts_rec_arrays; i++){
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
    double const Xe_current = Xe_ne_data.first;
    //check for non-zero and negative values
    double actXe_current = Xe_current < global_tol ? global_tol: Xe_current;
    Xe_arr_saha[i] = actXe_current;
  }

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

    // Initial Xe with be the last (valid) added value from the Saha regime. 
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

  }
  else{
    std::cout << "------ Xe-array completely filled by Saha equation. Are you sure this is correct? ------" << std::endl;
  }

  // Create new arrays containing the log of Xe and ne
  Vector log_Xe_arr(npts_rec_arrays);
  Vector log_XeSaha_arr(npts_rec_arrays);
  Vector log_ne_arr(npts_rec_arrays);

  for(int i=0; i<npts_rec_arrays; i++){
    log_Xe_arr[i] = log(Xe_arr[i]);
    log_XeSaha_arr[i] = log(Xe_arr_saha[i]);
    log_ne_arr[i] = log(ne_arr[i]);
  }

  // Make splines of the result
  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "Xe");
  log_XeSaha_of_x_spline.create(x_array, log_XeSaha_arr, "XeSaha");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");

  Utils::EndTiming("Xe");
}

std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  
  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  

  // Compute Xe from the Saha equation
  double Tb = cosmo->get_TCMB(x);
  double nb = get_nb_of_x(x);

  double b = 1./nb * std::pow(k_b*m_e*Tb/(2.0*M_PI*hbar*hbar), 1.5) * exp(-epsilon_0/(k_b*Tb));
  if(b>1e7)
    Xe = 1.0;
  else{
    Xe = (-b+sqrt(b*b+4.0*b))/2.0;
  }

  const double nH = get_nH_of_x(x);
  ne = Xe*nH;

  return std::pair<double,double>(Xe, ne);
}

int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Constants for finding RHS of peebles eq.
  const double Tb = cosmo->get_TCMB(x);
  const double eps_tb = epsilon_0/(k_b*Tb);
  const double alpha = c/hbar*sqrt(3.0*sigma_T/(8.0*M_PI))*m_e; // Dimensionless fine structure constant. 
  const double H = cosmo->H_of_x(x);

  // Finding the terms involved in the RHS
  const double phi2 = 0.448 * log(eps_tb); // dimensionless
  const double alpha2 = hbar*hbar/c*64.0*M_PI*alpha*alpha*phi2 / (m_e * m_e) * sqrt(eps_tb/(27.0*M_PI)); // dimension m^3/s
  const double beta = alpha2* std::pow(k_b*m_e*Tb/(2.0*M_PI*hbar*hbar), 1.5) * exp(-eps_tb);
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

  dXedx[0] = rhs;

  return GSL_SUCCESS;
}

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("Optical depth");

  const int npts = int(1e6);

  // Make reversed array that we integrate over
  Vector x_array_tau_rev = Utils::linspace(0, -x_start, npts);

  // The ODE system dtau/dx
  ODESolver tauOde;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Since we are using positive increasing x-values, we must provide the negative x-values in the functions for ne_of_x and H_of_x. Notice that this is actually -dtaudx in order for ut to integrate backwards. 
    dtaudx[0] = c*sigma_T*ne_of_x(-x)/(cosmo->H_of_x(-x));

    return GSL_SUCCESS;
  };

  // Initial Xe will be the last added value from the Saha regime. 
  Vector tau_init_vec{0.0};

  tauOde.solve(dtaudx, x_array_tau_rev, tau_init_vec);

  auto tau_vec_inv = tauOde.get_data_by_component(0);


  // Reverse the arrays
  Vector x_array_tau(npts);
  Vector tau_vec(npts);
  for(int i=0;i<npts;i++){
    x_array_tau[i] = -x_array_tau_rev[npts-1-i];
    tau_vec[i] = tau_vec_inv[npts-1-i];
  }

  // Generate dtaudx
  Vector dtau_vec(npts);
  for(int i=0;i<npts;i++){
    dtau_vec[i] = -c*ne_of_x(x_array_tau[i])*sigma_T / (cosmo->H_of_x(x_array_tau[i]));
  }

  // spline the results
  tau_of_x_spline.create(x_array_tau, tau_vec, "tau");
  dtaudx_of_x_spline.create(x_array_tau, dtau_vec, "dtaudx");

  Utils::EndTiming("Optical depth");
  Utils::StartTiming("Visibility function");

  // Generate gtilde
  Vector g_tilde(npts);
  for(int i=0;i<npts;i++){
    double const x_loc = x_array_tau[i];
    g_tilde[i] = -dtaudx_of_x(x_loc) * exp(-tau_of_x(x_loc));
  }


  // Generate dgtildedx
  Vector dg_tildedx(npts);
  for(int i=0;i<npts;i++){
    double x = x_array_tau[i];
    dg_tildedx[i] = exp(-tau_of_x(x)) * (dtaudx_of_x(x)*dtaudx_of_x(x) - ddtauddx_of_x(x));
  }
  
  g_tilde_of_x_spline.create(x_array_tau, g_tilde, "g");
  dg_tildedx_of_x_spline.create(x_array_tau, dg_tildedx, "dgdx");

  Utils::EndTiming("Visibility function");
}

void RecombinationHistory::solve_for_sound_horizon(){

  Utils::StartTiming("Sound horizon");

  const int npts = int(1e6);
  Vector x_array_s = Utils::linspace(x_start, x_end, npts);
  ODESolver sOde;
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    dsdx[0] = get_cs_of_x(x) / cosmo->Hp_of_x(x);
    return GSL_SUCCESS;
  };
  double init_val = get_cs_of_x(x_array_s[0])/cosmo->Hp_of_x(x_array_s[0]);
  Vector s_init_vec{init_val};
  sOde.solve(dsdx, x_array_s, s_init_vec);
  auto s_vec = sOde.get_data_by_component(0);
  s_of_x_spline.create(x_array_s, s_vec, "s");

  Utils::EndTiming("Sound horizon");

}

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return dtaudx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return dg_tildedx_of_x_spline(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return dg_tildedx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::XeSaha_of_x(double x) const{
  return exp(log_XeSaha_of_x_spline(x));
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

double RecombinationHistory::get_R_of_x(double x) const{
  return (3.0 * cosmo->get_OmegaB(x))/(4.0*cosmo->get_OmegaR(x));
}

double RecombinationHistory::get_cs_of_x(double x) const{
  double R = get_R_of_x(x);
  return c*sqrt(1./(3.*(1.+R)));
}

double RecombinationHistory::get_s_of_x(double x) const{
  return s_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "==================================================\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "==================================================\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = int(1e5);
  const double x_min   = -12.0;
  const double x_max   = 0;

  fp << "  x  , " << "   Xe  , " << " XeSaha  , " << "  ne   , " << " tau  , " << "  dtaudx    , " << " ddtauddx  , " << "  g    , " << "    dgdx  , " <<  "  ddgddx   ," << "\n";

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " , ";
    fp << Xe_of_x(x)           << " , ";
    fp << XeSaha_of_x(x)       << " , ";
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

void RecombinationHistory::analysis_output(const std::string filename) const{
  const int npts = int(1e6);
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  //  Find when last scattering happens (peak of visibility function)
  double gmax, x_LS;
  x_LS = x_array[0];
  gmax = g_tilde_of_x(x_LS);

  // Iterate through and save values for tau close to 1
  for(int i=0;i<npts; i++){
    double x_temp = x_array[i];
    double gtemp = g_tilde_of_x(x_temp);
    if(gtemp > gmax){
      x_LS = x_temp;
      gmax = g_tilde_of_x(x_LS);
    }
  }

  // Find out when recombination happened: Xe=0.1
  double Xe_min, XeSaha_min, x_rec, x_recSaha;
  x_rec = x_recSaha = x_array[0];
  Xe_min = Xe_of_x(x_rec);
  XeSaha_min = XeSaha_of_x(x_recSaha);

  for(int i=0; i<npts; i++){
    double x_temp = x_array[i];
    double Xe_temp = Xe_of_x(x_temp);
    double XeSaha_temp = XeSaha_of_x(x_temp);

    if(abs(Xe_temp-0.1) < abs(Xe_min-0.1)){
      x_rec = x_temp;
      Xe_min = Xe_of_x(x_rec);
    }
    if(abs(XeSaha_temp-0.1) < abs(XeSaha_min-0.1)){
      x_recSaha = x_temp;
      XeSaha_min = XeSaha_of_x(x_recSaha);
    }
  }
  // std::cout << "tau min: "<<tau_min<< "    "<<"x_LS: "<< x_LS<<std::endl;
  // std::cout << "Xe min: "<<Xe_min<< "    "<<"x_rec: "<< x_rec<<std::endl;
  // std::cout << "XeSaha min: "<<XeSaha_min<< "    "<<"x_recSaha: "<< x_recSaha<<std::endl;

  //  Get last scattering times
  double t_LS = cosmo->get_t_of_x(x_LS);
  double z_LS = cosmo->get_z_of_x(x_LS);

  //  Get recombination times
  double t_rec = cosmo->get_t_of_x(x_rec);
  double z_rec = cosmo->get_z_of_x(x_rec);

  //  Get Saha recombination times
  double t_recSaha = cosmo->get_t_of_x(x_recSaha);
  double z_recSaha = cosmo->get_z_of_x(x_recSaha);

  //  Get the sound horizons
  double rs_LS = get_s_of_x(x_LS);
  double rs_rec = get_s_of_x(x_rec);
  double rs_recSaha = get_s_of_x(x_recSaha);

  double freeze_out = Xe_of_x(0);
  double freeze_outSaha = XeSaha_of_x(0);

  std::cout <<"\n";
  std::cout << "==================================================\n";
  std::cout << "Analysis of recombination class: \n";
  std::cout << "==================================================\n";
  std::cout << "Phenomenon      " << "    x    " << "    z     "  << "     t  [Myr]" << "     r_s [Mpc] " << "\n";
  std::cout << "Last scattering " << x_LS << "  " << z_LS << "  " << t_LS*toMyr <<  "  " << rs_LS*toMpc << "\n";
  std::cout << "Recombination   " << x_rec << "  " << z_rec <<"  " << t_rec*toMyr <<  "  " << rs_rec*toMpc << "\n";
  std::cout << "Saha recomb     " << x_recSaha << "  " << z_recSaha <<"  " << t_recSaha*toMyr <<  "  " << rs_recSaha*toMpc << "\n";
  std::cout << "Freeze out      " << freeze_out <<"\n";
  std::cout << "Freeze out Saha " << freeze_outSaha << "\n";
  std::cout << std::endl;


  // write to file
  std::ofstream fp(filename.c_str());
  fp << "Phenomenon  " << ",";
  fp << "    x       " << ",";
  fp << "    z       " << ",";
  fp << "    t [Myr] " << ",";
  fp << "   r_s [Mpc]" << ",";
  fp << "\n";
  fp << "Last scattering " << " , " << x_LS << " , " << z_LS << " , " << t_LS*toMyr <<  " , " << rs_LS*toMpc <<" , " << "\n";
  fp << "Recombination   " << " , " << x_rec << " , " << z_rec << " , " << t_rec*toMyr <<  " , " << rs_rec*toMpc <<" , " << "\n";
  fp << "Saha            " << " , " << x_recSaha << " , " << z_recSaha << " , " << t_recSaha*toMyr << " , " << rs_recSaha*toMpc << " , " <<"\n";

}


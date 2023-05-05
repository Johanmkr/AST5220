#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

/*
  DO ALL THE SOLVING
*/

void PowerSpectrum::solve(){ 

  /*
    CREATE THE k-ARRAY WITH DIFFERENT RESOLUTION
  */

  double log_k_min = log(k_min);
  double log_k_max = log(k_max);

  // k-s for LOS integral
  double delta_k_LOS = 2.*M_PI / (eta0 * n_k_theta);
  Vector k_theta_array = get_linspace_from_delta(k_min, k_max, delta_k_LOS); // Linearly spaced logarithmic values
  Vector log_k_theta_array = log(k_theta_array);

  // k-s for power spectrum integral
  double delta_k_ps = 2.*M_PI / (eta0 * n_k_ps);
  Vector k_ps_array = get_linspace_from_delta(k_min, k_max, delta_k_ps);
  Vector log_k_ps_array = log(k_ps_array);


  /*
    GENERATE SPLINES OF BESSE FUNCTION
  */

  generate_bessel_function_splines();


  /*
    PERFORM THE LINE OF SIGHT INTEGRAL AND GENERATE SPLINE
  */

  line_of_sight_integration(k_theta_array);


  /*
    CALCULATE THE CMB POWER SPECTRUM AND SPLINE THE RESULT
  */

  Utils::StartTiming("cell_TT");
  auto cell_TT = solve_for_cell(log_k_ps_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  Utils::EndTiming("cell_TT");
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}


// Small utility functions

Vector PowerSpectrum::get_linspace_from_delta(double min, double max, double delta){
  double number_points = (int)abs(((max-min)/delta));
  return Utils::linspace(min, max, number_points);
}

double PowerSpectrum::get_finite_integral(Vector x_arr, Vector y_arr){
  // Declare and define needed variables 
  double integral_value = 0;
  for(size_t i=0; i<x_arr.size()-1; i++){
    integral_value += 0.5 * (y_arr[i]+y_arr[i+1]) * (x_arr[i+1]-x_arr[i]);
  }
  return integral_value;
}

/*
  Generate splines of j_ell(z) needed for LOS integration
*/

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");

  int const size_ell = (int)ells.size();
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(size_ell);
    
  // Determine argument interval
  double z_min    = 0.0;
  double z_max    = k_max*eta0;
  double delta_z  = 2.*M_PI/n_bessel;
  Vector z_array  = get_linspace_from_delta(z_min, z_max, delta_z);

  #pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < size_ell; i++){
    const int ell  = ells[i];

    //  Vector to store j_ell values
    Vector j_ell_arr(z_array.size());

    //  Loop across z-array
    for(size_t i=0; i<z_array.size(); i++){
      j_ell_arr[i] = Utils::j_ell(ell, z_array[i]);
    }
    j_ell_splines[i].create(z_array, j_ell_arr);
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  //  Find step length in x direction
  double delta_x = 2.*M_PI / n_LOS;

  // Create arrays
  Vector x_array = get_linspace_from_delta(x_start, x_end, delta_x);
  Vector integrand(x_array.size());
  
  for(size_t ik = 0; ik < k_array.size(); ik++){
    double const k_val = k_array[ik]; // k-value for each iteration
    for(size_t il = 0; il < ells.size(); il++){
      double const ell = ells[il]; // ell-value for each iteration

      // Quite ineffective with severeal foor-loops
      for(size_t i=0; i<x_array.size(); i++){
        integrand[i] = source_function(x_array[i], k_val) * j_ell_splines[il](k_val*(eta0 - cosmo->eta_of_x(x_array[i])));
      }
      result[il][ik] = get_finite_integral(x_array, integrand);
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  for(size_t il=0; il<ells.size(); il++){
    thetaT_ell_of_k_spline[il].create(k_array, thetaT_ell_of_k[il]);
  }

}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  Vector result(nells);

  //  Loop over and integrate for all ells
  for(int il = 0; il<nells; il++){
    double ell = ells[il];
    Vector integrand(log_k_array.size());
    double k_val;
    for(size_t i=0; i<log_k_array.size(); i++){
      k_val = exp(log_k_array[i]);
      integrand[i] = primordial_power_spectrum(k_val) * f_ell_spline[il](k_val)*g_ell_spline[il](k_val);
    }
    result[il] = 4.*M_PI * get_finite_integral(log_k_array, integrand);
  }
  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}

double PowerSpectrum::get_thetaT_ell_of_k_spline(const int il_idx, const double k) const{
  return thetaT_ell_of_k_spline[il_idx](k);
}





//====================================================
// Output the cells to file
//====================================================


void PowerSpectrum::output_theta(std::string filename) const{
  // Output of theta l for l={6, 100, 200, 500, 1000}
  std::ofstream fp(filename.c_str());
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, 10000); // Linearly spaced logarithmic values
  Vector k_array = exp(log_k_array);
  double c = Constants.c;
  double H0 = cosmo->H0;

  fp << "ckH , " << "T_6  , " << "T_100 , "<<"T_200 , "<< "T_500 , "<<"T_1000 , "<<"\n";
  auto print_data = [&] (const double k){
    fp << c*k/H0 << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[0], k) << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[1], k) << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[2], k) << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[3], k) << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[4], k) << " , ";
    // fp << thetaT_ell_of_k_spline[test_ell_idx[0]](k) << " , ";
    // fp << thetaT_ell_of_k_spline[test_ell_idx[1]](k) << " , ";
    // fp << thetaT_ell_of_k_spline[test_ell_idx[2]](k) << " , ";
    // fp << thetaT_ell_of_k_spline[test_ell_idx[3]](k) << " , ";
    // fp << thetaT_ell_of_k_spline[test_ell_idx[4]](k) << " , ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  fp << " ell , " << " cell_TT , " << "\n";
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " , ";
    fp << cell_TT_spline( ell ) * normfactor  << " , ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}


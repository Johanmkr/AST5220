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


  
  if(C_l_separation){
    std::cout <<"\n"<<"Performing LOS integration for all effects separately..."<<"\n"<<std::endl;
    //  For SW effect
    Utils::StartTiming("C_L separation");
    std::function<double(double,double)> source_function_T_SW = [&](double x, double k){
    return pert->get_Source_T_SW(x,k);
    };
    std::function<double(double,double)> source_function_T_ISW = [&](double x, double k){
    return pert->get_Source_T_ISW(x,k);
    };
    std::function<double(double,double)> source_function_T_DOP = [&](double x, double k){
    return pert->get_Source_T_DOP(x,k);
    };
    std::function<double(double,double)> source_function_T_POL = [&](double x, double k){
    return pert->get_Source_T_POL(x,k);
    };

    // SW
    line_of_sight_integration(k_theta_array, source_function_T_SW);
    auto cell_SW = solve_for_cell(log_k_ps_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
    Cell_SW_spline.create(ells, cell_SW, "Cell_SW_of_ell");

    // ISW
    line_of_sight_integration(k_theta_array, source_function_T_ISW);
    auto cell_ISW = solve_for_cell(log_k_ps_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
    Cell_ISW_spline.create(ells, cell_ISW, "Cell_ISW_of_ell");

    // DOP
    line_of_sight_integration(k_theta_array, source_function_T_DOP);
    auto cell_DOP = solve_for_cell(log_k_ps_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
    Cell_DOP_spline.create(ells, cell_DOP, "Cell_DOP_of_ell");

    // POL
    line_of_sight_integration(k_theta_array, source_function_T_POL);
    auto cell_POL = solve_for_cell(log_k_ps_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
    Cell_POL_spline.create(ells, cell_POL, "Cell_POL_of_ell");

    Utils::EndTiming("C_L separation");
  }

 // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };
  line_of_sight_integration(k_theta_array, source_function_T);

  /*
    CALCULATE THE CMB POWER SPECTRUM AND SPLINE THE RESULT
  */

  Utils::StartTiming("cell_TT");
  auto cell_TT = solve_for_cell(log_k_ps_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  Utils::EndTiming("cell_TT");
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
  Vector x_array = get_linspace_from_delta(-8.7, x_end, delta_x);
  
  // for(size_t ik = 0; ik < k_array.size(); ik++){
  //   double const k_val = k_array[ik]; // k-value for each iteration
  //   for(size_t il = 0; il < ells.size(); il++){
  //     double const ell = ells[il]; // ell-value for each iteration

  //     // Quite ineffective with severeal foor-loops
  //     for(size_t i=0; i<x_array.size(); i++){
  //       integrand[i] = source_function(x_array[i], k_val) * j_ell_splines[il](k_val*(eta0 - cosmo->eta_of_x(x_array[i])));
  //     }
  //     result[il][ik] = get_finite_integral(x_array, integrand);
  //   }
  // }


  // double k_val;
  // double ell;
  // size_t il;
  // size_t i;
  std::cout<<"Performing LOS integral..."<<std::endl;
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t ik = 0; ik < k_array.size(); ik++){
    if( (10*ik) /k_array.size() != (10*ik+10) / k_array.size()) {
        std::cout << (100*ik+100)/k_array.size() << "% " << std::flush;
        if(ik == k_array.size()-1) std::cout << std::endl;
      }
    double k_val = k_array[ik]; // k-value for each iteration
    for(size_t il = 0; il < ells.size(); il++){
      double ell = ells[il]; // ell-value for each iteration

      // double integrand_current = source_function(x_array[0], k_val) * j_ell_splines[il](k_val*(eta0-cosmo->eta_of_x(x_array[0])));
      // double integrand_next;

      // double integral_value = 0.0;

      // for(size_t i=0; i<x_array.size()-1; i++){
      //   integrand_next = source_function(x_array[i+1], k_val) * j_ell_splines[il](k_val*(eta0 - cosmo->eta_of_x(x_array[i+1])));
      //   integral_value += 0.5 * (integrand_current+integrand_next) * (x_array[i+1]-x_array[i]);
      //   integrand_current = integrand_next;
      // }
      // result[il][ik] = integral_value;

      Vector integrand(x_array.size());
      // Quite ineffective with several for-loops
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
void PowerSpectrum::line_of_sight_integration(Vector & k_array, std::function<double(double,double)> &source_function){
  const int n_k        = k_array.size();
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function);

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
// P(k) in units of
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k) const{
  // Variables and constants
  double c      = Constants.c;
  double Phi    = pert->get_Phi(x, k);
  double OmegaM = cosmo->get_OmegaM(0.0);
  double H0     = cosmo->H0;
  double a      = exp(x);
  double P_prim = primordial_power_spectrum(k) * 2. * M_PI * M_PI / ( k * k * k);

  // Calculate Delta_M
  double Delta_M = (2*c*c*k*k*Phi*a) / (3*OmegaM*H0*H0);

  // Calculate P(k,x)
  double pofk = Delta_M*Delta_M * P_prim ;
  
  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}  

double PowerSpectrum::get_cell_SW(const double ell) const{
  return Cell_SW_spline(ell);
}
double PowerSpectrum::get_cell_ISW(const double ell) const{
  return Cell_ISW_spline(ell);
}
double PowerSpectrum::get_cell_DOP(const double ell) const{
  return Cell_DOP_spline(ell);
}
double PowerSpectrum::get_cell_POL(const double ell) const{
  return Cell_POL_spline(ell);
}

double PowerSpectrum::get_thetaT_ell_of_k_spline(const int il_idx, const double k) const{
  return thetaT_ell_of_k_spline[il_idx](k);
}

double PowerSpectrum::get_source_func(const double x, const double k) const{
  return pert->get_Source_T(x, k);
}

double PowerSpectrum::get_bessel_func(const int il_idx, const double z) const{
  return j_ell_splines[il_idx](z);
}

double PowerSpectrum::get_LOS_integrand(const double x, const double k, const int il_idx) const{
  double bessel_argument = k * (eta0 - cosmo->eta_of_x(x));
  return get_source_func(x,k) * j_ell_splines[il_idx](bessel_argument);
}





/*
  WRITE OUTPUT TO CORRESPONDING FILE
*/

void PowerSpectrum::output_bessel(std::string filename) const{
  std::ofstream fp(filename.c_str());
  Vector z_array = Utils::linspace(0.0,1000,5000);
  fp << "z , " << "j_6 , "<<"j_100 , "<< "j_200 , " << "j_500 , " << "j_1000 , " << "\n";
  auto print_data = [&] (const double z) {
    fp << z << " , ";
    fp << get_bessel_func(test_ell_idx[0], z) << " , ";
    fp << get_bessel_func(test_ell_idx[1], z) << " , ";
    fp << get_bessel_func(test_ell_idx[2], z) << " , ";
    fp << get_bessel_func(test_ell_idx[3], z) << " , ";
    fp << get_bessel_func(test_ell_idx[4], z) << " , ";
    fp << "\n";
  };
  std::for_each(z_array.begin(), z_array.end(), print_data);
}

void PowerSpectrum::output_LOS_integrand(std::string filename) const{
  std::ofstream fp(filename.c_str());
  Vector x_array = Utils::linspace(x_start, x_end, 10000);
  fp << "x , " << "Sj_6k_min , " << " Sj_100k_min , " << "Sj_200k_min , " << "Sj_ 500k_min , " << "Sj_1000k_min , " << 
                  "Sj_6k_max , " << " Sj_100k_max , " << "Sj_200k_max , " << "Sj_ 500k_max , " << "Sj_1000k_max , " << "\n";
  auto print_data = [&] (const double x){
    fp << x << " , ";
    fp << get_LOS_integrand(x, k_min, test_ell_idx[0]) << " , ";
    fp << get_LOS_integrand(x, k_min, test_ell_idx[1]) << " , ";
    fp << get_LOS_integrand(x, k_min, test_ell_idx[2]) << " , ";
    fp << get_LOS_integrand(x, k_min, test_ell_idx[3]) << " , ";
    fp << get_LOS_integrand(x, k_min, test_ell_idx[4]) << " , ";
    fp << get_LOS_integrand(x, k_max, test_ell_idx[0]) << " , ";
    fp << get_LOS_integrand(x, k_max, test_ell_idx[1]) << " , ";
    fp << get_LOS_integrand(x, k_max, test_ell_idx[2]) << " , ";
    fp << get_LOS_integrand(x, k_max, test_ell_idx[3]) << " , ";
    fp << get_LOS_integrand(x, k_max, test_ell_idx[4]) << " , ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}


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
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

void PowerSpectrum::output_Cl_integrand(std::string filename) const{
  std::ofstream fp(filename.c_str());
  Vector k_array = Utils::linspace(k_min, k_max, 10000);
  
  fp << " k  , " << "  T_lK , " << "\n";
  auto print_data = [&] (const double k){
    fp << k << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[0], k) * get_thetaT_ell_of_k_spline(test_ell_idx[0], k) / k << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[1], k) * get_thetaT_ell_of_k_spline(test_ell_idx[1], k) / k << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[2], k) * get_thetaT_ell_of_k_spline(test_ell_idx[2], k) / k << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[3], k) * get_thetaT_ell_of_k_spline(test_ell_idx[3], k) / k << " , ";
    fp << get_thetaT_ell_of_k_spline(test_ell_idx[4], k) * get_thetaT_ell_of_k_spline(test_ell_idx[4], k) / k << " , ";
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

void PowerSpectrum::output_C_l_sep(std::string filename) const{
  if(C_l_separation){
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  fp << " ell , " << " cell_SW , " << " cell_ISW , "<< " cell_DOP , " << "cell_POL , " << "\n";
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    // double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
    //   * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    // double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " , ";
    fp << get_cell_SW(ell) * normfactor  << " , ";
    fp << get_cell_ISW(ell) * normfactor  << " , ";
    fp << get_cell_DOP(ell) * normfactor  << " , ";
    fp << get_cell_POL(ell) * normfactor  << " , ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
  }
  else{
    std::cout << "No separation" << std::endl;
  }
}


void PowerSpectrum::output_MPS(std::string filename) const{
  std::ofstream fp(filename.c_str());
  double Mpc = Constants.Mpc;
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, 10000);
  Vector k_array = exp(log_k_array);
  double h = cosmo->get_h();
  fp << "k ," << " Pk , " << "\n";
  auto print_data = [&] (const double k){
    fp << k * Mpc/ h << " , ";
    fp << get_matter_power_spectrum(0.0, k) * (h * h * h) / (Mpc*Mpc*Mpc) << " , ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}


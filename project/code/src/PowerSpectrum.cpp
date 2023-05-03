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

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){ 

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, n_k); // Linearly spaced logarithmic values
  Vector k_array = exp(log_k_array);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================

  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  Utils::StartTiming("cell_TT");
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
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

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");

  int const size_ell = (int)ells.size();
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(size_ell);
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================
  
  // Find z_max
  double z_max           = k_max * eta0 > 35000 ? 35000 : k_max * eta0; //not sure if this works

  //  Determine the number of points in the array based on 10 samples per oscillation
  double delta_z        = 2.0*M_PI / (50.);

  #pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < size_ell; i++){
    const int ell  = ells[i];

    // The range vary with l
    // double delta_zl   = delta_z / ell;
    int nr_z          = (int)(z_max / delta_z);
    // std::cout<<i<<" / "<<size_ell<<" - nr_z: "<<nr_z<<std::endl;
    // int nr_z          = 250000;
    Vector z_array    = Utils::linspace(0.0, z_max, nr_z);

    //  Vector to store j_ell values
    Vector j_ell_arr(nr_z);

    //  Loop across z-array
    for(int i=0; i<nr_z; i++){
      j_ell_arr[i] = Utils::j_ell(ell, z_array[i]);
      // j_ell_arr[i] = std::sph_bessel(ell, z_array[i]);
      // if(std::isnan(j_ell_arr[i])){
      //   std::cout<<ell<<std::endl;
      // }
      // std::cout<<j_ell_arr[i]<<std::endl;
    }
    // Make the j_ell_splines[i] spline
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
  
  for(size_t ik = 0; ik < k_array.size(); ik++){
    double const k_val = k_array[ik]; // k-value for each iteration
    double const delta_variabel = 2.*M_PI/(10.*k_val);  // Constant part of delta_x interval
    double const const_bessel_argument = k_val*eta0;  // constant part of bessel argument
    for(size_t il = 0; il < ells.size(); il++){
      double const ell = ells[il];

      //  Variables needed to perform the integral across x for the source function
      double x_current = x_start;
      double integral_sum = 0;
      double bessel_argument = const_bessel_argument - k_val*cosmo->eta_of_x(x_current);
      double integrand_current = source_function(x_current, k_val) * j_ell_splines[il](bessel_argument);

      //  Variables to fill
      double x_next;
      double delta_x=0;
      double integrand_next;

      while(x_current+delta_x<x_end){
        // delta_x = delta_variabel * cosmo->Hp_of_x(x_current);
        delta_x = 0.001;
        // std::cout << delta_x << std::endl;
        x_next = x_current+delta_x;
        bessel_argument = const_bessel_argument - k_val * cosmo->eta_of_x(x_next);
        integrand_next = source_function(x_next, k_val) * j_ell_splines[il](bessel_argument);

        // if(std::isnan(j_ell_splines[il](bessel_argument))){
        // std::cout<<"il: "<<il<<" - j_l: "<<j_ell_splines[il](bessel_argument)<<" - arg: "<<bessel_argument<<std::endl;
        // }

        integral_sum += 0.5 * (integrand_current + integrand_next) * delta_x;

        // Iterative scheme 
        x_current = x_next;
        integrand_current = integrand_next;
      }

      result[il][ik] = integral_sum;
      // std::cout<<integral_sum<<std::endl;
    }
    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    // ...
    // ...
    // ...

    // Store the result for Source_ell(k) in results[ell][ik]
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
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



  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  // ...
  // ...
  // ...
  // ...

  Vector result(nells);
  
  //  Loop over and integrate for all ells
  for(int il = 0; il<nells; il++){
    double ell = ells[il];
      // Declare and define variables for the integration
      double integral_sum = 0;
      double k_current = exp(log_k_array[0]);
      double integrand_current = primordial_power_spectrum(k_current) * f_ell_spline[il](k_current) * g_ell_spline[il](k_current);
      double integrand_next;
      double k_next;
      for(size_t i=0; i<log_k_array.size()-1; i++){
        k_next = exp(log_k_array[i+1]);
        // std::cout << k_next << std::endl;
        integrand_next = primordial_power_spectrum(k_next) * f_ell_spline[il](k_next) * g_ell_spline[il](k_next);
        // std::cout << f_ell_spline[il](k_next) << std::endl;
        integral_sum += 0.5 * (integrand_current + integrand_next) * (k_next - k_current);

        k_current = k_next;
        integrand_current = integrand_next;
      }
      result[il] = integral_sum;
      // std::cout << integral_sum << std::endl;
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


//====================================================
// Output the cells to file
//====================================================

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


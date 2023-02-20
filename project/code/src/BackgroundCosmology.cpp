#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //  H0
  H0 = h * Constants.H0_over_h; // FIXME: is this correct?

  // OmegaR
  // OmegaR = 5.047e-5 * std::pow(TCMB/2.7255, 4) * (0.7/h)*(0.7/h);
  OmegaR = 16 * M_PI * M_PI * M_PI * Constants.G * std::pow(Constants.k_b*TCMB, 4) / (90 * std::pow(Constants.hbar*Constants.c, 3)*Constants.c*Constants.c * H0 * H0);

  // OmegaNu
  // OmegaNu = 3.491e-5 * std::pow(TCMB/2.7255, 4) * (0.7/h) * (0.7/h) * (Neff/3.046);
  OmegaNu = Neff * 7./8. * std::pow(4./11., 4./3.)*OmegaR;

  //  OmegaLambda
  OmegaLambda = 1 - (OmegaB + OmegaCDM + OmegaK + OmegaR + OmegaNu);


  //  Total densities
  OmegaM = OmegaB+OmegaCDM;
  OmegaR_tot = OmegaNu + OmegaR;
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(int nr_points){
  // int nr_points = (int)1e4; // Change this accordingly

  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, nr_points);


  /*
  Making the eta spline:
  */

  // Utils::StartTiming("Eta");
  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
     detadx[0] = Constants.c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  Vector eta_initial{Constants.c/Hp_of_x(Constants.x_start)};
  ODESolver etax;
  etax.solve(detadx, x_array, eta_initial);
  Vector eta_array = etax.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "eta_of_x");

  // Utils::EndTiming("Eta");

  /*
  Making the t spline:
  */

  // Utils::StartTiming("t");
  // ODE for dtdx
  ODEFunction dtdx = [&](double x, const double *eta, double *dtdx){
    dtdx[0] = 1/H_of_x(x);
    return GSL_SUCCESS;
  };

  Vector t_initial{1/(2*H_of_x(Constants.x_start))};
  ODESolver tx;
  tx.solve(dtdx, x_array, t_initial);
  Vector t_array = tx.get_data_by_component(0);

  t_of_x_spline.create(x_array, t_array, "t_of_x");
  t0 = t_of_x_spline(0.0);
  // Utils::EndTiming("t");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  double Hx = H0 * sqrt(OmegaLambda + OmegaK*exp(-2*x) + (OmegaB + OmegaCDM)*exp(-3*x) + (OmegaR+OmegaNu)*exp(-4*x));
  return Hx;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  return H0 * sqrt(u_of_x(x));
}

double BackgroundCosmology::u_of_x(double x) const{
  return OmegaLambda*exp(2*x) + OmegaK + (OmegaB+OmegaCDM)*exp(-x) + (OmegaR+OmegaNu)*exp(-2*x);
}

double BackgroundCosmology::dudx_of_x(double x) const{
  return 2*OmegaLambda*exp(2*x) -(OmegaB+OmegaCDM)*exp(-x) -2*(OmegaR+OmegaNu)*exp(-2*x);
}

double BackgroundCosmology::dduddx_of_x(double x) const{
  return 4*OmegaLambda*exp(2*x) + (OmegaB+OmegaCDM)*exp(-x) +4*(OmegaR+OmegaNu)*exp(-2*x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  return H0/(2*sqrt(u_of_x(x))) * dudx_of_x(x);
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  return H0*(-1/4*std::pow(u_of_x(x), -3./4.) * dudx_of_x(x)*dudx_of_x(x) + 1/(2*sqrt(u_of_x(x)))*dduddx_of_x(x));
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  // XXX: Optimisable
  return OmegaB * H0 * H0 * exp(-3*x) / (H_of_x(x) * H_of_x(x));
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  // XXX: Optimisable
  return OmegaR * H0 * H0 * exp(-4*x) / (H_of_x(x) * H_of_x(x));
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;
  //  XXX: Optimisable
  return OmegaNu * H0 * H0 * exp(-4*x) / (H_of_x(x) * H_of_x(x));
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  //  XXX: Optimisable
  return OmegaCDM * H0 * H0 * exp(-3*x) / (H_of_x(x) * H_of_x(x));
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  //  XXX: Optimisable
  return OmegaLambda * H0 * H0 / (H_of_x(x) * H_of_x(x));
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  //  XXX: Optimisable
  return OmegaK * H0 * H0 * exp(-2*x) / (H_of_x(x) * H_of_x(x));
}

// EXTRA FUNCTIONs

double BackgroundCosmology::get_r_of_x(double x) const{
  // double OK = get_OmegaK(x);
  if(OmegaK==0.0) return get_comoving_distance_of_x(x);
  double r, argument = sqrt(abs(OmegaK))*H0*get_comoving_distance_of_x(x)/Constants.c;
  if(OmegaK<0){
    return sin(argument)/argument;
  }
  else{
    return sinh(argument)/argument;
  }
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  return get_r_of_x(x) * exp(x);
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  return get_r_of_x(x) * exp(-x);
  // return get_comoving_distance_of_x(x) * exp(-x);
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  // FIXME: check this.
  return Constants.c/Hp_of_x(Constants.x_start)-eta_of_x(x);
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "==================================================\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "==================================================\n";
  std::cout << "OmegaB:       " << OmegaB      << " "     <<"\n";
  std::cout << "OmegaCDM:     " << OmegaCDM    << " "     <<"\n";
  std::cout << "OmegaM:       " << OmegaM      << " "     <<"\n";
  std::cout << "OmegaLambda:  " << OmegaLambda << " "     <<"\n";
  std::cout << "OmegaK:       " << OmegaK      << " "     <<"\n";
  std::cout << "OmegaNu:      " << OmegaNu     << " "     <<"\n";
  std::cout << "OmegaR:       " << OmegaR      << " "     <<"\n";
  std::cout << "OmegaR_tot:   " << OmegaR_tot  << " "     <<"\n";
  std::cout << "Neff:         " << Neff        << " "     <<"\n";
  std::cout << "h:            " << h           << " "     <<"\n";
  std::cout << "TCMB:         " << TCMB        << " [K]"   <<"\n";
  std::cout << "t0:           " << t0/(1e9*365*24*60*60)<< " [Gyr]"  <<"\n";
  std::cout << "H0:           " << H0 << " [1/s]" <<"\n";
  std::cout << std::endl;
} 


void BackgroundCosmology::test() const{

  std::cout << "Sum of densities:         " << OmegaM + OmegaR_tot + OmegaLambda + OmegaK << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -20.0;
  const double x_max =  5.0;
  const int    n_pts =  1000;

  Vector x_array = Utils::linspace(x_min, x_max, n_pts);
  
  // Vector x1_array = Utils::linspace(-12.0, 0.0, n_pts);
  // Vector x2_array = Utils::linspace(-14.0, 0.0, n_pts);
  // Vector x3_array = Utils::linspace(-20.0, 5.0, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << get_luminosity_distance_of_x(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}


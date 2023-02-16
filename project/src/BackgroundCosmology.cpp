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

  // OmegaR
  // OmegaR = 5.047e-5 * std::pow(TCMB/2.7255, 4) * (0.7/h)*(0.7/h);
  OmegaR = 2 * M_PI * M_PI / 30 * std::pow(Constants.k_b*TCMB, 4) / (std::pow(Constants.hbar*Constants.c, 3)*Constants.c*Constants.c);

  // OmegaNu
  // OmegaNu = 3.491e-5 * std::pow(TCMB/2.7255, 4) * (0.7/h) * (0.7/h) * (Neff/3.046);
  OmegaNu = Neff * 7./8. * std::pow(4./11., 4./3.)*OmegaR;

  //  OmegaLambda
  OmegaLambda = 1 - (OmegaB + OmegaCDM + OmegaK + OmegaR + OmegaNu);

  //  H0
  H0 = h * Constants.H0_over_h; // FIXME: is this correct?
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  int nr_points = (int)1e4; // Change this accordingly

  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, nr_points);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
     detadx[0] = Constants.c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  Vector initial{0.0};
  ODESolver etax;
  etax.solve(detadx, x_array, initial);
  auto eta = etax.get_data();
  Utils::EndTiming("Eta");
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
  //  XXX: Optimisable
  return OmegaK * H0 * H0 * exp(-2*x) / (H_of_x(x) * H_of_x(x));
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  return get_comoving_distance_of_x(x) * exp(-x);
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  // chi = eta0-eta (what is eta0)?

  return 0.0;
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
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

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
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}


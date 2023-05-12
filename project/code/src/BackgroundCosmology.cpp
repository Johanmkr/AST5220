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
  H0 = h * Constants.H0_over_h;

  // OmegaR
  double k_bTCMB = Constants.k_b * TCMB;
  OmegaR = 16. * M_PI * M_PI * M_PI * Constants.G * (k_bTCMB*k_bTCMB*k_bTCMB*k_bTCMB)/ (90. * (Constants.hbar*Constants.hbar*Constants.hbar) * (Constants.c*Constants.c*Constants.c*Constants.c*Constants.c) * H0 * H0);

  // OmegaNu
  OmegaNu = Neff * 7./8. * std::pow(4./11., 4./3.)*OmegaR;

  //  OmegaLambda
  OmegaLambda = 1. - (OmegaB + OmegaCDM + OmegaK + OmegaR + OmegaNu);


  //  Total densities
  OmegaM = OmegaB+OmegaCDM;
  OmegaR_tot = OmegaNu + OmegaR;

  // Equalities
  x_RM = std::log(OmegaR_tot/OmegaM);
  x_ML = std::log(OmegaM/OmegaLambda)/3.;
  x_acc = std::log(OmegaM/(2.*OmegaLambda))/3.;

}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(int nr_points){
  Vector x_array = Utils::linspace(x_start, x_end+5, nr_points);

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
    dtdx[0] = 1./H_of_x(x);
    return GSL_SUCCESS;
  };

  Vector t_initial{1./(2.*H_of_x(Constants.x_start))};
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
  double Hx = H0 * sqrt(OmegaLambda + OmegaK*exp(-2.*x) + (OmegaB + OmegaCDM)*exp(-3.*x) + (OmegaR+OmegaNu)*exp(-4.*x));
  return Hx;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  return H0 * sqrt(u_of_x(x));
}

double BackgroundCosmology::u_of_x(double x) const{
  return OmegaLambda*exp(2.*x) + OmegaK + (OmegaB+OmegaCDM)*exp(-x) + (OmegaR+OmegaNu)*exp(-2.*x);
}

double BackgroundCosmology::dudx_of_x(double x) const{
  return 2.*OmegaLambda*exp(2.*x) -(OmegaB+OmegaCDM)*exp(-x) -2.*(OmegaR+OmegaNu)*exp(-2.*x);
}

double BackgroundCosmology::dduddx_of_x(double x) const{
  return 4.*OmegaLambda*exp(2.*x) + (OmegaB+OmegaCDM)*exp(-x) +4.*(OmegaR+OmegaNu)*exp(-2.*x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  return H0/(2.*sqrt(u_of_x(x))) * dudx_of_x(x);
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  return H0*(-1./(4.*(sqrt(u_of_x(x))*sqrt(u_of_x(x))*sqrt(u_of_x(x))))* dudx_of_x(x)*dudx_of_x(x) + 1./(2.*sqrt(u_of_x(x)))*dduddx_of_x(x));
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  double Hret = H0 / Hp_of_x(x);
  double alpha = 1.0;
  return OmegaB * exp(-alpha * x) * Hret * Hret;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  double Hret = H0 / Hp_of_x(x);
  double alpha = 2.0;
  return OmegaR * exp(-alpha * x) * Hret * Hret;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;
  double Hret = H0 / Hp_of_x(x);
  double alpha = 2.0;
  return OmegaNu * exp(-alpha * x) * Hret * Hret;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  double Hret = H0 / Hp_of_x(x);
  double alpha = 1.0;
  return OmegaCDM * exp(-alpha * x) * Hret * Hret;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  double Hret = H0 / Hp_of_x(x);
  double alpha = -2.0;
  return OmegaLambda * exp(-alpha * x) * Hret * Hret;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;    
  double Hret = H0 / Hp_of_x(x);
  double alpha = 0.0;
  return OmegaK * exp(-alpha * x) * Hret * Hret;
}

double BackgroundCosmology::get_OmegaM(double x) const{
  return get_OmegaB(x) + get_OmegaCDM(x);
}

double BackgroundCosmology::get_OmegaRtot(double x) const{
  return get_OmegaR(x) + get_OmegaNu(x);
}

// EXTRA FUNCTIONs

double BackgroundCosmology::get_r_of_x(double x) const{
  // double OK = get_OmegaK(x);
  if(abs(OmegaK) < 1e-5) return get_comoving_distance_of_x(x);
  double argument = sqrt(abs(OmegaK))*H0*get_comoving_distance_of_x(x)/Constants.c;
  if(OmegaK<0){
    return get_comoving_distance_of_x(x) * sin(argument)/argument;
  }
  else{
    return get_comoving_distance_of_x(x) * sinh(argument)/argument;
  }
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  return get_r_of_x(x) * exp(x);
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  return get_r_of_x(x) * exp(-x);
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  return eta_of_x(0)-eta_of_x(x);
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::deta_of_x(double x) const{
  return eta_of_x_spline.deriv_x(x);
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

double BackgroundCosmology::get_t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_z_of_x(double x) const{
  return exp(-x) - 1;
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
  std::cout << "t0:           " << t0*toGyr << " [Gyr]"  <<"\n";
  std::cout << "H0:           " << H0/toGyr << " [1/Gyr]" <<"\n";
  std::cout << std::endl;
} 


//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename, double x_min, double x_max, int n_pts) const{
  // const double x_min = -20.0;
  // const double x_max =  5.0;
  // const int    n_pts =  1000;

  Vector x_array = Utils::linspace(x_min, x_max, n_pts);
  
  // Vector x1_array = Utils::linspace(-12.0, 0.0, n_pts);
  // Vector x2_array = Utils::linspace(-14.0, 0.0, n_pts);
  // Vector x3_array = Utils::linspace(-20.0, 5.0, n_pts);

  std::ofstream fp(filename.c_str());
  fp << "   x       "                  << " , ";
  fp << "    eta     "        << " , ";
  fp << "    deta     "        << " , ";
  fp << "    Hp      "        << " , ";
  fp << "    dHp     "        << " , ";
  fp << "    ddHp    "    << " , ";
  fp << "    OmegaB  "      << " , ";
  fp << "    OmegaCDM"    << " , ";
  fp << " OmegaLambda" << " , ";
  fp << "    OmegaR  "      << " , ";
  fp << "    OmegaNu "     << " , ";
  fp << "    OmegaK  "      << " , ";
  fp << "    d_L     " << " , ";
  fp << "    t     " << " , ";
  fp <<"\n";
  auto print_data = [&] (const double x) {
    fp << x                  << " , ";
    fp << eta_of_x(x)        << " , ";
    fp << deta_of_x(x)        << " , ";
    fp << Hp_of_x(x)         << " , ";
    fp << dHpdx_of_x(x)      << " , ";
    fp << ddHpddx_of_x(x)    << " , ";
    fp << get_OmegaB(x)      << " , ";
    fp << get_OmegaCDM(x)    << " , ";
    fp << get_OmegaLambda(x) << " , ";
    fp << get_OmegaR(x)      << " , ";
    fp << get_OmegaNu(x)     << " , ";
    fp << get_OmegaK(x)      << " , ";
    fp << get_luminosity_distance_of_x(x) << " , ";
    fp << get_t_of_x(x) << " , ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

void BackgroundCosmology::write_table_of_important_values(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  fp << "Quantity" << " , " << "x" << " , " << "z" << " , " << "t [Gyr]" << " , " << "\n";
  fp << "RM-equality" << " , " << x_RM << " , " << get_z_of_x(x_RM) << " , " << get_t_of_x(x_RM)*toGyr << " , " << "\n";
  fp << "ML-equality" << " , " << x_ML << " , " << get_z_of_x(x_ML) << " , " << get_t_of_x(x_ML)*toGyr << " , " << "\n";
  fp << "Accel. start" << " , " << x_acc << " , " << get_z_of_x(x_acc) << " , " << get_t_of_x(x_acc)*toGyr << " , " << " \n";
  fp << "Age of universe" << " , " << x0 << " , " << get_z_of_x(x0) << " , " << get_t_of_x(x0)*toGyr << " , " << "\n";
  fp << "Conformal time" << " , " << x0 << " , " << get_z_of_x(x0) << " , " << eta_of_x(x0)/Constants.c*toGyr << " , " << "\n";
}


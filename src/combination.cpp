#include "combination.h"

using namespace comb;

Combination::Combination(double a0, double a12, double b0, double b12, double smax, double step_size)
:
a0{a0},
a12{a12},
b0{b0},
b12{b12},
smax{smax},
step_size{step_size}
{readin(); spline();}

double Combination::kaellen(double a, double b, double c)
{
	return std::pow(a,2)+std::pow(b,2)+std::pow(c,2)-2.*(a*b+b*c+a*c);
}

double Combination::t(double s, double z)
{ 
	double delta = std::pow(constants::mass_kaon(),2)*(std::pow(constants::mass_kaon(),2)-std::pow(constants::mass_pi(),2));
	double lambda = std::sqrt(kaellen(s,0,std::pow(constants::mass_kaon(),2)))*std::sqrt(kaellen(s,std::pow(constants::mass_pi(),2),std::pow(constants::mass_kaon(),2)));
	return 1./2*(3*s0-s+(lambda*z-delta)/s);
}

double Combination::u(double s, double z)
{
	return 3*s0 - s - t(s,z);
}

void Combination::readin()
{
	std::string directory = "../../gammaKKpi_amp/basisfunctions/";
	std::vector<std::string> v = {"f0_a.txt", "f0_b.txt", "f0_c.txt", "f12_a.txt", "f12_b.txt", "f12_c.txt", "G0.txt", "Gp.txt"};

	array_s.resize(8);
	array_real.resize(8);
	array_imag.resize(8);

    int i = 0;
    for(auto &file: v){
    	std::ifstream myfile(directory + file);
    	double s, abs, arg;
	    while (myfile >> s >> abs >> arg) {
	    	Complex comp = abs*std::exp(I*arg);
	        array_s[i].push_back(s);
		    array_real[i].push_back(std::real(comp));
		    array_imag[i].push_back(std::imag(comp));
	    }
	    i+=1;
    }
    return;
}

void Combination::spline()
{
	F0_a_real = gsl::Interpolate(array_s[0],array_real[0],gsl::InterpolationMethod::cubic);
	F0_b_real = gsl::Interpolate(array_s[1],array_real[1],gsl::InterpolationMethod::cubic);
	F0_c_real = gsl::Interpolate(array_s[2],array_real[2],gsl::InterpolationMethod::cubic);
	F12_a_real = gsl::Interpolate(array_s[3],array_real[3],gsl::InterpolationMethod::cubic);
	F12_b_real = gsl::Interpolate(array_s[4],array_real[4],gsl::InterpolationMethod::cubic);
	F12_c_real = gsl::Interpolate(array_s[5],array_real[5],gsl::InterpolationMethod::cubic);
	G0_real = gsl::Interpolate(array_s[6],array_real[6],gsl::InterpolationMethod::cubic);
	Gp_real = gsl::Interpolate(array_s[7],array_real[7],gsl::InterpolationMethod::cubic);
	F0_a_imag = gsl::Interpolate(array_s[0],array_imag[0],gsl::InterpolationMethod::cubic);
	F0_b_imag = gsl::Interpolate(array_s[1],array_imag[1],gsl::InterpolationMethod::cubic);
	F0_c_imag = gsl::Interpolate(array_s[2],array_imag[2],gsl::InterpolationMethod::cubic);
	F12_a_imag = gsl::Interpolate(array_s[3],array_imag[3],gsl::InterpolationMethod::cubic);
	F12_b_imag = gsl::Interpolate(array_s[4],array_imag[4],gsl::InterpolationMethod::cubic);
	F12_c_imag = gsl::Interpolate(array_s[5],array_imag[5],gsl::InterpolationMethod::cubic);
	G0_imag = gsl::Interpolate(array_s[6],array_imag[6],gsl::InterpolationMethod::cubic);
	Gp_imag = gsl::Interpolate(array_s[7],array_imag[7],gsl::InterpolationMethod::cubic);
	return;
}

Complex Combination::F0(double s)
{
	return a0*(F0_a_real(s)+I*F0_a_imag(s))+b0*(F0_b_real(s)+I*F0_b_imag(s))+F0_c_real(s)+I*F0_c_imag(s);
}

Complex Combination::F12(double s)
{
	return a12*(F12_a_real(s)+I*F12_a_imag(s))+b12*(F12_b_real(s)+I*F12_b_imag(s))+F12_c_real(s)+I*F12_c_imag(s);
}

Complex Combination::Gp(double s)
{
	return Gp_real(s) + I*Gp_imag(s);
}

Complex Combination::G0(double s)
{
	return G0_real(s) + I*G0_imag(s);
}

Complex Combination::f(int i, double s)
{
	Complex F, Fhat;
	if (i == 1){
		F = F12(s) - F0(s);
		auto integrand = [this, s](double z){return 3./4*(1-std::pow(z,2))*(Gp(t(s,z)) - G0(t(s,z)) + F12(u(s,z)) - F0(u(s,z)));};
		Fhat = cauchy::complex_integration(integrand, -1, 1);
	}
	else if (i == 2){
		F = -std::sqrt(2)*(F12(s) - F0(s));
		auto integrand = [this, s](double z){return 3./4*(1-std::pow(z,2))*std::sqrt(2)*(G0(t(s,z)) + F12(u(s,z)) + F0(u(s,z)));};
		Fhat = cauchy::complex_integration(integrand, -1, 1);
	}
	else if (i == 3){
		F = F12(s) + F0(s);
		auto integrand = [this, s](double z){return 3./4*(1-std::pow(z,2))*(Gp(t(s,z)) + G0(t(s,z)) + F12(u(s,z)) + F0(u(s,z)));};
		Fhat = cauchy::complex_integration(integrand, -1, 1);
	}
	else if (i == 4){
		F = std::sqrt(2)*(F12(s) + F0(s));
		auto integrand = [this, s](double z){return 3./4*(1-std::pow(z,2))*std::sqrt(2)*(G0(t(s,z)) - F12(u(s,z)) + F0(u(s,z)));};
		Fhat = cauchy::complex_integration(integrand, -1, 1);
	}
	else{
		throw std::domain_error("Function not defined for i!={1,2,3,4}.");
	}
	
	return F + Fhat;
}

void Combination::output(int i)
{
	std::string directory = "../../gammaKKpi_amp/";
	std::vector<std::string> v {"F1.dat", "F2.dat", "F3.dat", "F4.dat"};

	std::ofstream myfile(directory + v[i-1]);
	if (myfile.is_open())
  	{
  		for(double s=std::pow(constants::mass_kaon()+constants::mass_pi(),2); s < smax; s+=step_size)
  		{
  			myfile << s << "\t" << std::real(f(i,s)) << "\t" << std::imag(f(i,s)) << std::endl;
  		}	
	   	myfile.close();
  	}

	return;
}
#include "input.h"

using namespace input;

gammaKKpi::gammaKKpi(int i, bool use_err)
:
matchpoint{std::pow(1.,2)}, //in GeV^2
combination{0.9,1.0,-0.4,2.7,2,0.001} //Subtraktionconstants are set here
{which_readin(i, use_err); spline(); spline_err(); match();}

void gammaKKpi::which_readin(int i, bool use_err)
{
	if(use_err){
		readin_old(i);
	}
	else{
		readin(i);
	}
}

void gammaKKpi::readin(int i)
{
	std::string file;
	if(i==1){
		file = "../../gammaKKpi_amp/F1.dat";
	}
	else if(i==2){
		file = "../../gammaKKpi_amp/F2.dat";
	}
	else if(i==3){
		file = "../../gammaKKpi_amp/F3.dat";
	}
	else if(i==4){
		file = "../../gammaKKpi_amp/F4.dat";
	}
	else{
		throw std::domain_error("Function not defined for i!={1,2,3,4}.");
	}
	
    double s, real, imag;
    std::ifstream testfile(file);
    if(!testfile){
    	combination.output(i);
    }
    std::ifstream myfile(file);
    while (myfile >> s >> real >> imag) {
        s_list.push_back(s);
        real_list.push_back(real);
        imag_list.push_back(imag);
        real_err_list.push_back(real); //err not included
        imag_err_list.push_back(imag); //err not included
    }
    return;
}

void gammaKKpi::readin_old(int i)
{
	std::string file;
	if(i==1){
		file = "../../gammaKKpi_amp/F1_2sub.txt";
	}
	else if(i==2){
		file = "../../gammaKKpi_amp/F2_2sub.txt";
	}
	else if(i==3){
		file = "../../gammaKKpi_amp/F3_2sub.txt";
	}
	else if(i==4){
		file = "../../gammaKKpi_amp/F4_2sub.txt";
	}
	else{
		throw std::domain_error("Function not defined for i!={1,2,3,4}.");
	}

    double s, real, imag, real_err, imag_err, cross, cross_err;
    std::ifstream myfile(file);
    while (myfile >> s >> real >> imag >> real_err >> imag_err >> cross >> cross_err) {
        s_list.push_back(s);
        real_list.push_back(real);
        imag_list.push_back(imag);
        real_err_list.push_back(real_err);
        imag_err_list.push_back(imag_err);
    }
    return;
}

void gammaKKpi::spline()
{
	for(std::size_t i=0; i<real_list.size(); i++){
		abs_list.push_back(std::sqrt(std::pow(real_list[i],2)+ std::pow(imag_list[i],2)));
	}
	spline_abs = gsl::Interpolate(s_list,abs_list,gsl::InterpolationMethod::cubic);
}


void gammaKKpi::spline_err()
{
	for(std::size_t i=0; i<real_list.size(); i++){
		abs_err_list.push_back(std::sqrt(std::pow(real_list[i]*real_err_list[i]/abs_list[i],2)+ std::pow(imag_list[i]*imag_err_list[i]/abs_list[i],2)));
	}
	spline_abs_err = gsl::Interpolate(s_list,abs_err_list,gsl::InterpolationMethod::cubic);
}


void gammaKKpi::match()
{
	double step = 10e-3;
	double y = spline_abs(matchpoint);
	double y_err = spline_abs_err(matchpoint);
	std::pair<double, double> m;
	std::pair<double, double> m_err;
	m = gsl::derivative(spline_abs, matchpoint, step);
	m_err = gsl::derivative(spline_abs_err, matchpoint, step);
	
	a = -pow(y,2)/m.first;
	b = -matchpoint-y/m.first;
	a_err = -pow(y_err,2)/m_err.first;
	b_err = -matchpoint-y_err/m_err.first;
}


double gammaKKpi::cont(double s)
{
	return a/(s+b); //we use this function to continue with 1/s this leads to the definition of the parameters a and b in match()
}

double gammaKKpi::cont_err(double s)
{
	return a_err/(s+b_err); //we use this function to continue with 1/s this leads to the definition of the parameters a and b in match()
}


double gammaKKpi::operator()(double s)
{
	if(s <= matchpoint)
	{
		return spline_abs(s);
	}
	else if(s > matchpoint)
	{
		return cont(s);
	}
	else
	{
		throw std::domain_error{"s value is out of range!"};
	}
}

double gammaKKpi::operator[](double s)
{
	if(s <= matchpoint)
	{
		return spline_abs_err(s);
	}
	else if(s > matchpoint)
	{
		return cont_err(s);
	}
	else
	{
		throw std::domain_error{"s value is out of range!"};
	}
}
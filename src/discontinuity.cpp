#include "discontinuity.h"

using namespace disc;

Discontinuity::Discontinuity(bool charged_kaon_int)
:
charged_kaon_int{charged_kaon_int},
sth{std::pow(constants::mass_pi()+constants::mass_kaon(),2.)},
absgammaKKpic{input::gammaKKpi(1,false)},
absgammaKKpin{input::gammaKKpi(2,false)}
{}

double Discontinuity::lambda(double s, double m1, double m2) 
{
	return std::pow(s,2) + std::pow(m1,4) + std::pow(m2,4) - 2. * (s * std::pow(m1,2) + s * std::pow(m2,2) + std::pow(m1*m2,2));
}

double Discontinuity::operator()(double s)
{
	if(s<sth)
	{
		return 0;
	}
	else
	{
		if (charged_kaon_int){
			return 1./(4.*constants::pi()*std::pow(constants::elementary_charge(),2)) * std::pow(lambda(s,constants::mass_pi(), constants::mass_kaon()), 3./2) 
					 /(72*std::pow(s,2.)) * std::pow(absgammaKKpic(s),2.);
			// don't include the i here
		}
		else{
			return 1./(4.*constants::pi()*std::pow(constants::elementary_charge(),2)) * std::pow(lambda(s,constants::mass_pi(), constants::mass_kaon()), 3./2) 
					 /(72*std::pow(s,2.)) * std::pow(absgammaKKpin(s),2.);
			// don't include the i here
		}

	}
}

double Discontinuity::operator[](double s)
{
	if(s<sth)
	{
		return 0;
	}
	else
	{
		if (charged_kaon_int){
			return 1./(4.*constants::pi()*std::pow(constants::elementary_charge(),2)) * std::pow(lambda(s,constants::mass_pi(), constants::mass_kaon()), 3./2) 
					 /(72*std::pow(s,2.)) * 2 *  absgammaKKpic(s) * absgammaKKpic[s];
			// don't include the i here
		}
		else{
			return 1./(4.*constants::pi()*std::pow(constants::elementary_charge(),2)) * std::pow(lambda(s,constants::mass_pi(), constants::mass_kaon()), 3./2) 
					 /(72*std::pow(s,2.)) * 2 *  absgammaKKpin(s) * absgammaKKpin[s];
			// don't include the i here
		}

	}
}
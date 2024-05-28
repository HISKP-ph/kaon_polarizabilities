#include "dispersiveintegral.h"

using namespace disp;

DispersiveIntegral::DispersiveIntegral(int num_sub, int err)
:
gammaKKpicdisc{disc::Discontinuity(true)},
gammaKKpindisc{disc::Discontinuity(false)},
sth{gammaKKpicdisc.sth},
multiply{1e3},
subtraction_point{std::pow(constants::mass_kaon(),2)},
num_sub{num_sub},
err{err}
{set_cutoff();}

void DispersiveIntegral::set_cutoff(){
	if (num_sub==1){
		cutoff = std::numeric_limits<double>::infinity();
	}
	else if (num_sub > 1){
		cutoff = std::numeric_limits<double>::infinity();
	}
	else{
		throw std::domain_error("num_sub must be an integer greater or equal to 1.");
	}
}

double DispersiveIntegral::numerator(double s){
	if (err == 0){
		return (gammaKKpicdisc(s) + gammaKKpindisc(s))*multiply; 
	}
	else if (err == 1){
		return (gammaKKpicdisc(s) - gammaKKpicdisc[s]  + gammaKKpindisc(s) - gammaKKpindisc[s])*multiply; 
	}
	else if (err == 2){
		return (gammaKKpicdisc(s) + gammaKKpicdisc[s]  + gammaKKpindisc(s) + gammaKKpindisc[s])*multiply; 
	}
	else{
		throw std::domain_error("err must be 0 (without error), 1 (low) or 2 (up).");
	}
}

double DispersiveIntegral::integrand_cauchy(double s, double s_prime){
       return (numerator(s_prime)-numerator(s))/(std::pow(s_prime-subtraction_point,num_sub)*(s_prime-s)); //number of subtractions
}

double DispersiveIntegral::numeric_integral_cauchy(double s, double lower_limit){
    double mandelstam_s=s;

    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_cauchy(mandelstam_s,s_prime); }};
    
    auto integration = gsl::Cquad();
    integration.reserve(10000);
    integration.set_absolute(0.0);
    integration.set_relative(1e-7);

    gsl::Value result_0 = integration(integrand, lower_limit, cutoff);
    return std::get<0>(result_0);
}

double DispersiveIntegral::integrand_trivial(double s, double s_prime){
       return numerator(s_prime)/(std::pow(s_prime-subtraction_point,num_sub)*(s_prime-s)); //number of subtractions
} 

double DispersiveIntegral::numeric_integral_trivial(double s, double lower_limit){
    double mandelstam_s=s;
    const auto integrand{[mandelstam_s,this](double s_prime){
        return this->integrand_trivial(mandelstam_s,s_prime); }};

    auto integration = gsl::Cquad();
    integration.reserve(10000);
    integration.set_absolute(0.0);
    integration.set_relative(1e-7);

    gsl::Value result = integration(integrand, lower_limit, cutoff);
    return std::get<0>(result);
}

double DispersiveIntegral::integral_analytic(double s, double lower_limit){
	if (subtraction_point > lower_limit){
		throw std::domain_error("Subtraction_point must be smaller than lower_limit!");
	}
	if (num_sub==1){
		if(cutoff == std::numeric_limits<double>::infinity() ) {
			return -numerator(s) * std::log((s-lower_limit)/(lower_limit - subtraction_point));
		}
		else{
			return numerator(s) * (std::log((cutoff-s)/(cutoff-subtraction_point)) - std::log((s-lower_limit)/(lower_limit - subtraction_point)) );
		}
		
	}
	else if (num_sub == 2){
		if (cutoff == std::numeric_limits<double>::infinity()){
			return numerator(s)* (-(s-subtraction_point)/(lower_limit-subtraction_point) - std::log((s-lower_limit)/(lower_limit - subtraction_point)));
		}
		else{
			throw std::domain_error("2 sutbractions for cutoff != infinty not implemented.");
		}
	}
	else{
		throw std::domain_error("Numer of subtraction not implemented.");
	}
}

Complex DispersiveIntegral::operator()(double s){
	if (s > cutoff)
	{
		throw std::domain_error("s larger than cutoff of integral.");
	}
	else if (s < sth){
		return 1./multiply*(std::pow(s-subtraction_point,num_sub)/(2*constants::pi())*numeric_integral_trivial(s, sth));
	}
	else{
		return 1./multiply*(std::pow(s-subtraction_point,num_sub)/(2*constants::pi())*numeric_integral_cauchy(s, sth) + 1./(2*constants::pi()) * integral_analytic(s, sth) + 1.i/2. * numerator(s) ) ;
	}

}
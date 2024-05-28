#ifndef _dispersionintegral_
#define _dispersionintegral_

#include <iostream>
#include <cmath>
#include <complex>
#include "gsl_interface.h"
#include "cauchy.h"
#include "discontinuity.h"

/// Namespace to calculate dispersive integral for gamma K to gamma K 
namespace disp{

using type_aliases::Complex;
using namespace std::complex_literals;

/// Class to calculate dispersive integral for gamma K to gamma K
class DispersiveIntegral
{
public:
	DispersiveIntegral(int num_sub, int err);
	///<@param num_sub number of subtractions used for the dispersion integral. Must be greater or equal to 1. Only 1 and 2 are implemented.
	///<@param err integer that determines which function for the disconinuity are used. 0: normal function is evaluated, 1: uncertainty below, 2: uncertainty above

	/// Discontinuity for charged kaon in intermediate state
	disc::Discontinuity gammaKKpicdisc;
	/// Discontinuity for neutral kaon in intermediate state
	disc::Discontinuity gammaKKpindisc;

	/// Threshold set in Constructor via DispersiveIntegral::gammaKKpicdisc
	double sth;
	/// Number set in Constructor to increase numerical stability
	double multiply;
	/// Cutoff of the disperison integral, set via DispersiveIntegral::set_cutoff according to num_sub
	double cutoff;
	/// Subtraction point of the dispersion integral, set in Constructor to M_K^2 
	double subtraction_point;
	int num_sub;
	int err;

	/// Function to set cutoff
	void set_cutoff();

	/// Constructs numerator for the disperion integral
	double numerator(double s);
	/// Integrand for the Cauchy integral
	double integrand_cauchy(double s, double s_prime);
	/// Calculates integral for the Cauchy integral
	double numeric_integral_cauchy(double s, double lower_limit);
	/// Integrand for the trivial integral
	double integrand_trivial(double s, double s_prime);
	/// Integrates trivially (only valid for s<sth)
	double numeric_integral_trivial(double s, double lower_limit);
	/// Analytic part of the integral when s>sth and using the Cauchy integration. Only implemented for num_sub={1,2}
	double integral_analytic(double s, double lower_limit);
	/// returns the disperion integral. Uses the right combination of integrals according to where s lies
	Complex operator()(double s);

};

}

#endif
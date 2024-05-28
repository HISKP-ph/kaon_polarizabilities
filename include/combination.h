#ifndef _combination_
#define _combination_

#include <fstream>
#include <iostream>
#include <vector>
#include "constants.h"
#include "type_aliases.h"
#include "gsl_interface.h"
#include "cauchy.h"

using type_aliases::Complex;

/// Namespace for the combination of basis functions for gamma K to K pi from https://inspirehep.net/literature/1835296
namespace comb{

/// Class for the combination of basis functions for gamma K to K pi from https://inspirehep.net/literature/1835296
class Combination
{
public:
	Combination(double a0, double a12, double b0, double b12, double smax, double step_size);
	///<@param a0 Subtraction constant for F^(0) amplitude
	///<@param a12 Subtraction constant for F^(1/2) amplitude
	///<@param b0 Subtraction constant for F^(0) amplitude
	///<@param b12 Subtraction constant for F^(1/2) amplitude
	///<@param smax Maximum value of Mandelstam s for the output file
	///<@param step_size Step size of s for the output file

	double s0 = 1./3*(2*std::pow(constants::mass_kaon(),2)+std::pow(constants::mass_pi(),2));
	double kaellen(double a, double b, double c);
	///< Källén function
	double t(double s,double z);
	///< Mandelstam t depending on
	///< @param s Mandelstam s
	///< @param z cosine of the s-channel scattering angle
	double u(double s,double z);
	///< Mandelstam u depending on
	///< @param s Mandelstam s
	///< @param z cosine of the s-channel scattering angle
	Complex I{0,1};

	Complex f(int i, double s);
	///< partial wave amplitude depending on
	///<@param s Mandelstam s
	///<@param i integer that maps to the notation in https://inspirehep.net/literature/1835296 as follows
	///< i=1: -0; i=2: 0-; i=3: 00; i=4: -+

	/// combines the basis functions to F^(0)
	Complex F0(double s);
	/// combines the basis functions to F^(1/2)
	Complex F12(double s);
	/// combines the basis functions to G^(+)
	Complex Gp(double s);
	/// combines the basis functions to G^(0)
	Complex G0(double s);

	double a0, a12, b0, b12;
	double smax, step_size;

	std::vector<std::vector<double>> array_s, array_real, array_imag;

	/// Called in Contructor. Reads the files for the basis functions; check the directory!
	void readin();
	/// Called in Contructor. Interpolates the basis functions with gsl 
	void spline();
	/// Called in Contructor. Outputs the partial wave corresponding to i (see definition of Combination::f ) into file 
	void output(int i);

	gsl::Interpolate F0_a_real;
	gsl::Interpolate F0_b_real;
	gsl::Interpolate F0_c_real;
	gsl::Interpolate F12_a_real;
	gsl::Interpolate F12_b_real;
	gsl::Interpolate F12_c_real;
	gsl::Interpolate Gp_real;
	gsl::Interpolate G0_real;
	gsl::Interpolate F0_a_imag;
	gsl::Interpolate F0_b_imag;
	gsl::Interpolate F0_c_imag;
	gsl::Interpolate F12_a_imag;
	gsl::Interpolate F12_b_imag;
	gsl::Interpolate F12_c_imag;
	gsl::Interpolate Gp_imag;
	gsl::Interpolate G0_imag;
};

}


#endif
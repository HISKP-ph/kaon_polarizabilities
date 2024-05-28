#ifndef _input_
#define _input_

#include <fstream>
#include <iostream>
#include <vector>
#include "gsl_interface.h"
#include "combination.h"

/// Namespace to readin, spline and match the gamma K to K pi function from https://inspirehep.net/literature/1835296
namespace input{

/// Class to readin, spline and match the gamma K to K pi function from https://inspirehep.net/literature/1835296
class gammaKKpi
{
public:
	gammaKKpi(int i, bool use_err);
	///<@param i integer that maps to the notation in https://inspirehep.net/literature/1835296 as follows
	///< i=1: -0; i=2: 0-; i=3: 00; i=4: -+
	///<@param use_err if true readin_old is used and files with uncertainties generated in https://inspirehep.net/literature/1835296 are used
	///< if true readin is used and new files using basis functions from https://inspirehep.net/literature/1835296 but adaptive subtraction constants are used, no uncertainties are calculated

	std::vector<double> s_list;
	std::vector<double> real_list;
	std::vector<double> imag_list;
	std::vector<double> real_err_list;
	std::vector<double> imag_err_list;
	std::vector<double> abs_list;
	std::vector<double> abs_err_list;
	double a;
	double b;
	double a_err;
	double b_err;
	/// Matchpoint set in Constructor
	double matchpoint;

	/// Called in gammaKKpi::readin if file does not exist with the subtraction constants calculated in https://inspirehep.net/literature/1835296
	comb::Combination combination;

	gsl::Interpolate spline_abs;
	gsl::Interpolate spline_abs_err;

	/// Called in Constructor. Depending on use_err decides which readin function is used
	void which_readin(int i, bool use_err);
	/// Reads in file. If the file does not exist creates new one according to gammaKKpi::combination
	void readin(int i);
	/// Reads in file from https://inspirehep.net/literature/1835296 including uncertainties
	void readin_old(int i);
	/// Splines functions
	void spline();
	/// Splines uncertainies
	void spline_err();

	/// Calculates the derivative. Needed for matching.
	double deriv(double x, double step);
	/// Matches the interpolated function to analytic function at gammaKKpi::matchpoint
	void match();

	/// Analytic function for continuation. Uses parameters from gammaKKpi::match
	double cont(double s);
	/// Analytic function for continuation of uncertainies. Uses parameters from gammaKKpi::match
	double cont_err(double s);

	/// Evalutes the function. Below the matchpoint uses interpolated splines, above the function gammaKKpi::cont
	double operator()(double s);
	/// Evalutes the uncertanty of the function. Below the matchpoint uses interpolated splines, above the function gammaKKpi::cont_err
	double operator[](double s);
};

}

#endif
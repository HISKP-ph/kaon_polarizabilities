#ifndef _discontinuity_
#define _discontinuity_

#include "constants.h"
#include "input.h"

/// Namespace to calculate the discontinuity for kaon polarizabilities
namespace disc{

/// Class to calculate the discontinuity for kaon polarizabilities
class Discontinuity
{
public:
	Discontinuity(bool charged_kaon_int);
	///<@param charged_kaon_int if true uses gamma K to K pi amplitude with charged kaon in the final state -> charged kaon in intermediate state for gamma K to gamma K
	///< if fales uses neutral kaon

	bool charged_kaon_int;

	/// Threshold for gamma K to K pi process
	double sth;

	/// gamma K to K pi amplitude with charged kaon in final state
	input::gammaKKpi absgammaKKpic;
	/// gamma K to K pi amplitude with neutral kaon in final state
	input::gammaKKpi absgammaKKpin;


	double lambda(double s, double m1, double m2);
	///< Källén function for masses
	///<@param s Mandelstam s
	///<@param m1 input for mass 1, not mass squared!
	///<@param m2 input for mass 2, not mass squared!

	/// Operator that evaluates the discontinuity
	double operator()(double s);
	/// Operator that evaluates the uncertainity of the discontinuity
	double operator[](double s);
};

}

#endif
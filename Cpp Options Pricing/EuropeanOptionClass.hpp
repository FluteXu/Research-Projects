// EuropeanOptionClass.hpp
// Functions and class implementation based on ClassforExerciseA

#ifndef CLASSFOREXERCISEA_HPP
#define CLASSFOREXERCISEA_HPP

#include <string>
#include <boost/math/distributions/normal.hpp>

#include "UtilitiesDJD/VectorsAndMatrices/Vector.cpp"

using namespace std;

// Global function to calculate call price (put price is much the same)
//double CallPrice(double T, double K, double sig, double r, double S);

class EuropeanOption
{
private:		

	void init();	// Initialise all default values
	void copy(const EuropeanOption& source);

	// 'Kernel' functions for option calculations
	double CallPrice() const;  //All parameter are encapsulated in member data of the class
	double PutPrice() const;
	double CallDelta() const;
	double PutDelta() const;


	// Overloading function for using divided differences method
    double CallPrice(double U) const;
	double PutPrice(double U) const;

	// Function using divided differences to calculate Delta
	double DiffCallDelta() const;
	double DiffPutDelta() const;
	

	// Gaussian functions
	double n(double x) const;
	double N(double x) const;


public:

	// Member data public for convenience; anyway, the format will 
	// not change for a plain option.

	double T;		// Expiry date
	double K;		// Strike price
	double sig;		// Volatility
	double r;		// Interest rate
	double b;       // Cost-of-carry
	double S;       // Current stock price

	string optType;	// Option name (call, put)


public:	// Public functions
	EuropeanOption();							// Default call option
	EuropeanOption(double T_input, double K_input, double sig_input,
		double r_input, double b_input, double S_input, string optType_input);  // Manual constructor
	EuropeanOption(const EuropeanOption& source);	// Copy constructor
	virtual ~EuropeanOption();	

	EuropeanOption& operator = (const EuropeanOption& source);

	// Functions that calculate option price and sensitivities
	double Price() const;
	double Delta() const;

	// Function using divided differences to calculate Delta
	double DiffDelta() const;


	//Selectors
	void Print() const;

	// Modifier functions
	void toggle();		// Change option type (C/P, P/C)


};

// Function produces a mesh array for A.1 d-e
// Function pointer to member function was used here
Vector<double, long> MeshArray(double Start, double End, double Interval, 
	string element_vary, EuropeanOption opt, double (EuropeanOption::*func)() const);


#endif
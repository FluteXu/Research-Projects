// EuropeanOptionClass.cpp
// Methods modified on ClassforExerciseA
//
// Modification dates:
// 2016-08-14 Christine Xu

#include "EuropeanOptionClass.hpp"
#include <cmath>
#include <iostream>

using namespace std;

//double CallPrice(double T, double K, double sig, double r, double S)
//{
//	using namespace boost::math;
//    normal_distribution<> StandardNormal(0, 1.0);
//
//	double tmp = sig * sqrt(T);
//
//	double d1 = ( log(S/K) + (r + (sig*sig)*0.5 ) * T )/ tmp;
//	double d2 = d1 - tmp;
//
//	return (S * cdf(StandardNormal, d1)) - (K * exp(-r * T)* cdf(StandardNormal, d2));
//}

//Implement normal distribution probability function using boost
double EuropeanOption::n(double x) const
{
	using namespace boost::math;
    normal_distribution<> StandardNormal(0, 1.0);

	return pdf(StandardNormal, x);
}

double EuropeanOption::N(double x) const
{
	using namespace boost::math;
    normal_distribution<> StandardNormal(0, 1.0);

	return cdf(StandardNormal, x);
}

// Kernel Functions (Haug)
double EuropeanOption::CallPrice() const
{
	double tmp = sig * sqrt(T);

	double d1 = ( log(S/K) + (b + (sig*sig)*0.5 ) * T )/ tmp;
	double d2 = d1 - tmp;

	return (S * exp((b-r)*T) * N(d1)) - (K * exp(-r * T)* N(d2));
}

double EuropeanOption::PutPrice() const
{
	double tmp = sig * sqrt(T);

	double d1 = ( log(S/K) + (b + (sig*sig)*0.5 ) * T )/ tmp;
	double d2 = d1 - tmp;

	return (K * exp(-r * T)* N(-d2)) - (S * exp((b-r)*T) * N(-d1));
}

// Overloading function for using divided differences method
double EuropeanOption::CallPrice(double U) const
{

	double tmp = sig * sqrt(T);

	double d1 = ( log(U/K) + (b+ (sig*sig)*0.5 ) * T )/ tmp;
	double d2 = d1 - tmp;


	return (U * exp((b-r)*T) * N(d1)) - (K * exp(-r * T)* N(d2));

}

double EuropeanOption::PutPrice(double U) const
{

	double tmp = sig * sqrt(T);
	double d1 = ( log(U/K) + (b+ (sig*sig)*0.5 ) * T )/ tmp;
	double d2 = d1 - tmp;

	return (K * exp(-r * T)* N(-d2)) - (U * exp((b-r)*T) * N(-d1));

}


// Functions calculate Delta
double EuropeanOption::CallDelta() const
{
	double tmp = sig * sqrt(T);

	double d1 = ( log(S/K) + (b + (sig*sig)*0.5 ) * T )/ tmp;


	return exp((b-r)*T) * N(d1);
}

double EuropeanOption::PutDelta() const
{
	double tmp = sig * sqrt(T);

	double d1 = ( log(S/K) + (b + (sig*sig)*0.5 ) * T )/ tmp;

	return exp((b-r)*T) * (N(d1) - 1.0);
}



// Function using divided differences to calculate Delta
double EuropeanOption::DiffCallDelta() const
{
	double h = 0.0001;
	return (CallPrice(S+h) - CallPrice(S-h)) / (2*h);
}
double EuropeanOption::DiffPutDelta() const
{
	double h = 0.0001;
	return (PutPrice(S+h) - PutPrice(S-h)) / (2*h);
}

/////////////////////////////////////////////////////////////////////////////////////

void EuropeanOption::init()
{	// Initialise all default values

	// Default values
	r = 0.05;
	b = 0;
	sig= 0.2;
	K = 110.0;
	T = 0.5;
	S = 100;
	optType = "C";		// European Call Option (this is the default type)
}

void EuropeanOption::copy( const EuropeanOption& source)
{
	r	= source.r;
	b   = source.b;
	sig = source.sig;	
	K	= source.K;
	T	= source.T;
	S	= source.S;
	
	optType = source.optType;
}

EuropeanOption::EuropeanOption() // Defualt constructor
{
	init();
}

EuropeanOption::EuropeanOption(double T_input, double K_input, double sig_input,
		double r_input, double b_input, double S_input, string optType_input)  // Manual constructor
{
	r	= r_input;
	b   = b_input;
	sig = sig_input;	
	K	= K_input;
	T	= T_input;
	S	= S_input;
	
	optType = optType_input;
}

EuropeanOption::EuropeanOption(const EuropeanOption& source)  // Conpy constructor
{
	copy(source);
}

EuropeanOption::~EuropeanOption()
{
}

EuropeanOption& EuropeanOption::operator = (const EuropeanOption& source)
{
	if (this == &source) return *this;

	copy (source);

	return *this;
}

double EuropeanOption::Price() const
{
	if (optType == "C")
	{	
		return CallPrice();
	}
	else
	{
		return PutPrice();
	}
}

double EuropeanOption::Delta() const
{
	if (optType == "C")
		return CallDelta();
	else
		return PutDelta();
}

// Function using divided differences to calculate Delta
double EuropeanOption::DiffDelta() const
{
	if (optType == "C")
		return DiffCallDelta();
	else
		return DiffPutDelta();
}

// Selectors
void EuropeanOption::Print() const
{
	cout << "T "<< T << ", K " << K << ", sig " << sig << ", r "
		<< r << ", b " << b <<", S " << S << endl;
}

// Modifier functions
void EuropeanOption::toggle()
{ // Change option type (C/P, P/C)

	if (optType == "C")
		optType = "P";
	else
		optType = "C";
}

// Function produces a mesh array for A.1 d-e
Vector<double, long> MeshArray(double Start, double End, double Interval, 
	string element_vary, EuropeanOption opt, double (EuropeanOption::*func)() const)
{
	int size;
	size = int((End+Interval-Start)/Interval);
	Vector<double, long> result(size, 0); 
	
	// Determine which elements would be vaired
	double* temp_pt;
	double orig;
	if (element_vary == "S")
	{
		temp_pt = &opt.S;
		orig = opt.S;
	}
	else if (element_vary == "T")
	{
		temp_pt = &opt.T;
		orig = opt.T;
	}
	else if (element_vary == "sig")
	{
		temp_pt = &opt.sig;
		orig = opt.sig;
	}


	// Make the corresponding element vary and store the price result
	for (long i = 0; i < size; i++)
	{
		*temp_pt = Start + i*Interval;
		result[i] = (opt.*func)();
	}

	*temp_pt = orig;  // Restore the original value
	

	return result;
}



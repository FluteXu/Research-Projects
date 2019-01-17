// TestOptionPrice.cpp
// Modified on TestSingleCurve.cpp to print the output in Excel
//
// Modification dates:
// 2016-08-14 Chrisntine Xu


#include "UtilitiesDJD/VectorsAndMatrices/Vector.cpp"
#include "UtilitiesDJD/ExcelDriver/ExcelMechanisms.hpp"
#include "UtilitiesDJD/Geometry/Range.cpp"
#include "UtilitiesDJD/ExceptionClasses/DatasimException.hpp"
#include "EuropeanOptionClass.hpp"

#include <cmath>
#include <list>
#include <string>

int main()
{
	// DON'T FORGET TO MODIFY EXCELIMPORTS.CPP for correct version of Excel.

	// Create EuropeanOption
	EuropeanOption Option = EuropeanOption(0.25, 65, 0.3, 0.08, 0.08, 60, "C");

	// Create PriceMeshArray to print
	Vector<double, long> PriceMeshArray = MeshArray(60, 70, 1, "S", Option, &EuropeanOption::Price);

	// Create S array
	double Start = 60; 
	double End = 70; 
	double Interval = 1;
	int size = int((End+Interval-Start)/Interval);
	Vector<double, long> S(size, 0);  // Default start index is 1

	for (long i = 0; i < size; i++)
	{
		S[i] = Start + i*Interval;
	}
	
	std::cout << "Data has been created\n";


	// Print the output to Excel
	try 
	{
		printOneExcel(S, PriceMeshArray, std::string("Option Price"), std::string("S"), std::string("Price"));
	}
	catch(DatasimException& e)
	{
		e.print();
	}

	//try 
	//{
	//	printOneExcel(S, PriceMeshArray); // Default annotation settings
	//}
	//catch(DatasimException& e)
	//{
	//	e.print();
	//}

	return 0;
}

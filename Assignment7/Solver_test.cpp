//
//#####################################################################
//                            Solver_test.cpp 
//#####################################################################
//
//
// This is a routine that uses a procedure of Non-pivotal Gaussian Elimination
// to solve systems of tri-diagonal equations.
// It accepts as input the coefficients of the tri-diagonal matrix in
// three vectors (one for the a's, one for the b's, and one for the c's)
// and the vector corresponding to the right hand side of the equation, z.
// The routine returns the vector x that is the solution to Px = z,
// where P is the tri-diagonal matrix.
//
// Output is printed to files of the form
//
//         uCompXXX.dat 
//
// where XXX is the output file index.
//
// Currently the format of the data is an ASCII file
// and can be read by either GNUplot or Matlab.
//
// Modify the program if you want Excel output.
//
//#####################################################################
// Created for Math 269B. 
// Version: Feb. 22, 2010
// Author : Yasser Taima
// From a template by : Chris Anderson
//#####################################################################
#include <iostream>
using namespace std;

#include <math.h>
#include <stdio.h>

#include "DoubleArray1D.h"         // 1D array class 

#include "ProfileOutput.h"         // Utility class for outputting data files in 
                                   // an appropriate format. 

void TridiagSolver(DoubleArray1D&, DoubleArray1D&, DoubleArray1D&, DoubleArray1D&, DoubleArray1D&);

int main(int argc, char* argv[]) 
{
    double xA       = 0.0;      // left endpoint
    double xB       = 1.0;      // right endpoint

    long N;                     	    // size of matrix

    long outputIndex;			    // index of output steps

    DoubleArray1D uN;                   // Array to hold computed solution values
    DoubleArray1D a;                        // Array to hold diagonal values
    DoubleArray1D b;                        // Array to hold lower diagonal values
    DoubleArray1D c;			    // Array to hold upper diagonal values
    DoubleArray1D z;			    // Array to hold input vector
    DoubleArray1D x;			    // Array to hold solution

    const char*  outputFileNamePrefix = "uComp";   // output filename prefix.
    int  outputFormat;                              // to specify output format


    cout << " Enter size of the matrix N : " << endl;
    cin  >> N;
 
    uN.initialize(N);  // N data points for the computed solution
    a.initialize(N);  // N data points in the matrix
    b.initialize(N);  // N-1 data points on the lower diagonal; b(0) is a dummy
    c.initialize(N); // N-1 data points on the lower diagonal; c(N-1) is a dummy
    z.initialize(N);
    x.initialize(N);
    //
    // Instantiate output class
    //

    ProfileOutput profileOut(xA,xB);

    // Test data, simplified
    long i;
    i = 0;
    // set a, b, c
    cout << " Enter diagonal element a : " << endl;
    cin  >> a(i);
    cout << " Enter lower diagonal element b : " << endl;
    cin  >> b(i);
    cout << " Enter upper diagonal element c : " << endl;
    cin  >> c(i);
    cout << " Enter z element : " << endl;
    cin  >> z(i);

    // set arrays a, b ,c
    for (i=1; i <= N-1; i++)
    {
	a(i) = a(i-1);
	b(i) = b(i-1);
	c(i) = c(i-1);
	z(i) = z(i-1)+1.0;
    }

    TridiagSolver(a, b, c, z, x);
    uN = x;

    //
    // Output the initial solution profile
    //
	// Output formats supported
	//
	// ProfileOutput::MATLAB
	// ProfileOutput::GNUPLOT
	// ProfileOutput::EXCEL
	// 
	// 
 

    outputIndex = 0;
    outputFormat = ProfileOutput::GNUPLOT;
    profileOut.output(uN, outputIndex, outputFileNamePrefix, outputFormat);
 
    return 0;
 }

//
// This routine accepts as input the coefficients of the tri-diagonal matrix
// in three vectors (one for the a's, one for the b's, and one for the c's)
// and the vector corresponding to the right hand side of the equation.
// The routine returns the vector x that is the solution to Px = z.

// Note that the b(0) and the c(N-1) elements are not used. The lengths of all
// vectors have been kept at N, for clarity of the code doing the computation.

void TridiagSolver(DoubleArray1D& aIn, DoubleArray1D& bIn, DoubleArray1D& cIn, DoubleArray1D& zIn, DoubleArray1D& xOut)
{
	DoubleArray1D e;
	DoubleArray1D d;
	DoubleArray1D y;
	DoubleArray1D x;

	long N = aIn.getSize(); // matrix size = size of the diagonal

	e.initialize(N);
	d.initialize(N);
	y.initialize(N);
	x.initialize(N);

	long i;
//
// Set e0 
//
	i = 0;
	e(i) = aIn(i);
	d(i) = 0.0;
//        printf(" e = %10.5f  Timestep : %ld \n", e(i), i);
//
// Set d and e
// Note that the b(0) and the c(N-1) elements are not used. The lengths of all
// vectors have been kept at N, for clarity of the code doing the computation.

	for(i = 1; i <= N-1; i++)
	{
		d(i) = bIn(i)/e(i-1);
//        printf("d = %10.5f  Timestep : %ld \n", d(i), i);
		e(i) = aIn(i) - d(i)*cIn(i-1);
//        printf("e = %10.5f  Timestep : %ld \n", e(i), i);
	}

// forward substitution
	i = 0;
	y(i) = zIn(i);

        for(i = 1; i <= N-1; i++)
        {
                y(i) = zIn(i) - d(i)*y(i-1);
//        printf("y = %10.5f  Timestep : %ld \n", y(i), i);
        }

// Backward substitution
        i = N-1;
	x(i) = y(i)/e(i);
//        printf("x = %10.5f  Timestep : %ld \n", x(i), i);

        for(i = N-2; i >= 0; i--)
        {

                x(i) = (y(i) - cIn(i)*x(i+1))/e(i);
//        printf("x = %10.5f  Timestep: %ld \n", x(i), i);
        }

//
	for(i = 0; i <= N-1; i++)
	{
		
		xOut(i) = x(i);
	}
}

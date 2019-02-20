//
//#####################################################################
//                            Assign7.cpp 
//#####################################################################
//
// This program computes the exact solution to 
// 
// u_t + a u_xx = 0 
//
// for xA <= x <= xB, 0 <= t <= tFinal with 
// u(0,x) = u0(x), u(0,t) = g0(t) and u(1,t) = g1(t) 
// boundary conditions. 
//
// The spatial mesh is uniform with M panels between
// x = xA and x = xB. The array indices run from 0 to M 
// as indicated below: 
// 
//
//               M panels 
// xA                                  xB
//   
// |---x---x---x---x---x---x---x---x---x
// 0   1   2   3   4   5               M  
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
// Version: Feb. 24, 2010
// Author : Chris Anderson 
// Modified by : Yasser Taima
//#####################################################################
#include <iostream>
using namespace std;

#include <math.h>
#include <stdio.h>

#include "DoubleArray1D.h"         // 1D array class 
#include "SinExactPeriodicSolution.h" // Class for the exact solution 

#include "ProfileOutput.h"         // Utility class for outputting data files in 
                                   // an appropriate format. 

void advance(DoubleArray1D&, DoubleArray1D&, double, double, double);
double maxNorm(DoubleArray1D&, double);
void TridiagSolver(DoubleArray1D&, DoubleArray1D&, DoubleArray1D&, DoubleArray1D&, DoubleArray1D&);

int main(int argc, char* argv[]) 
{
    double tFinal;              // final time

    double xA       = 0.0;      // left endpoint
    double xB       = 1.0;      // right endpoint
    double a        = -0.1;      // diffusion constant

    long M;                     // number of panels
    double dx;                  // mesh spacing 

    double    time;             // total simulation time
    double      dt;             // timestep
    double  dtTemp;             // timestep for sub-stepping 
    double  uErr_norm;

    long   outputCount;         // number of timsteps output
    double outputTimeIncrement; // time interval between output times
    double outputTime;          // next output time
    long   outputIndex;         // index of output steps

    DoubleArray1D uExact;       // Array to hold Exact solution values
    DoubleArray1D uN;           // Array to hold computed solution values
    DoubleArray1D uNp1;         // Array to hold computed next step solution values
    DoubleArray1D uErr;

    SinExactPeriodicSolution exactSoln(xA,xB,a);         // Exact solution class instance
    const char*  outputExactFileNamePrefix = "uExact";   // Exact output filename prefix.
    const char*  outputFileNamePrefix = "uComp";         // output filename prefix.
    int  outputFormat;                                   // to specify output format


    cout << " Enter number of grid panels M : " << endl;
    cin  >> M;
 
    cout <<  "Enter Final Time : " << endl;
    cin  >>  tFinal;

    cout <<  "Number of Output Times : " << endl;
    cin  >>  outputCount;

    //
    // Initialize spatial variables 
    //

    dx = (xB-xA)/double(M);

    dt = dx;

    uErr_norm = 0.0;

    //
    // Initialize time-stepping variables
    //

    time                  = 0.0;
    outputTimeIncrement   = tFinal/double(outputCount);
    outputTime            = time + outputTimeIncrement;
    outputIndex           = 0;
  
    //
    // Initialize solution using the exact solution evaluated at
    // t = tInitial 
    //

    uExact.initialize(M+1);  // M+1 data points in our domain with M panels
    uN.initialize(M+1);      // M+1 data points in our domain with M panels
    uNp1.initialize(M+1);    // M+1 data points in our domain with M panels
    uErr.initialize(M+1);
    exactSoln.evaluateExactSolution(uExact,0.0);
    exactSoln.evaluateExactSolution(uN,0.0);

    //
    // Instantiate output class
    //

    ProfileOutput profileOut(xA,xB);

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
    outputFormat = ProfileOutput::GNUPLOT;
    profileOut.output(uN, outputIndex, outputFileNamePrefix, outputFormat);
    profileOut.output(uExact, outputIndex, outputExactFileNamePrefix, outputFormat);
 
    printf("Output %3ld Time = %10.5f  Timesteps taken : %d \n", outputIndex, time,0);
    
    // Main time-stepping loop 

    int outputFlag       = 0;
    long kStep           = 0;
    int exitFlag	 = 0;			// XXXXX Correction 1/31/10
    while((time < tFinal)&&(exitFlag == 0))	// XXXXX Correction 1/31/10
    {

    // Simulation evolution 

    if(time + dt < outputTime)              // next step requires no output   
    {
        //
        // Advance numerical solution 
        //
	advance(uN, uNp1, dx, dt, a);
	uN = uNp1;

        // update time, stepcount, and exact solution 
        time += dt;
        kStep++;
        exactSoln.evaluateExactSolution(uExact,time); 
    }
    else                                    // next step requires output
    {
        outputFlag = 1;                      
        dtTemp     = outputTime - time;     // substep so we hit output time exactly

        if(dtTemp > 1.0e-10)                // advance the solution to output time
        {                                   // (if necessary) 
        //
        // Advance numerical solution 
        //
	advance(uN, uNp1, dx, dtTemp, a);
	uN = uNp1;
        
	// update time and exact solution 

        time += dtTemp;
        kStep++;
        exactSoln.evaluateExactSolution(uExact,time);
        }
	if(fabs(time-tFinal) < 1.0e-10){exitFlag = 1;} // XXXXX Correction 1/31/10

    }

    // Simulation output 

    if(outputFlag  == 1)
    {
        outputFlag  = 0;                    // reset flags and time markers
        outputTime += outputTimeIncrement;
        outputIndex++;      
        
        uErr = uN - uExact;
        uErr_norm = maxNorm(uErr, dx);

        printf("Output %3ld Time = %10.5f  Timesteps taken : %ld \n", outputIndex, time, kStep);

        outputFormat = ProfileOutput::GNUPLOT;
        profileOut.output(uN, outputIndex, outputFileNamePrefix, outputFormat);
        profileOut.output(uExact, outputIndex, outputExactFileNamePrefix, outputFormat);
        printf("The max norm of the error vector of length %3ld at time %5.2f is %10.5f \n", M, time, uErr_norm);
 
    }

    }

    return 0;
 }

//
// This routine advances the solution to u_t + a u_xx = 0 one timestep using
// a method based upon a trapezoidal rule time derivative approximation and
// a two-sided difference spatial derivative approximation.
// 
//
void advance(DoubleArray1D& uIn, DoubleArray1D& uOut, double h, double k, double a)
{
	long M = uIn.getSize() -1; // Number of panels = number of points - 1
	long i;

//  Arrays for the tridiagonal solver
	DoubleArray1D a_diag;                   // Array to hold diagonal values
	DoubleArray1D b;                        // Array to hold lower diagonal values
	DoubleArray1D c;                        // Array to hold upper diagonal values
	DoubleArray1D z;                        // Array to hold input vector
	DoubleArray1D x;                        // Array to hold solution
	
	a_diag.initialize(M-1);  // N data points in the matrix
	b.initialize(M-1);  // N-1 data points on the lower diagonal; b(0) is a dummy
	c.initialize(M-1);  // N-1 data points on the lower diagonal; c(N-1) is a dummy
	z.initialize(M-1); 
	x.initialize(M-1); 
//
//  Left edge point --- use boundary condition
//
	i = 0;
	uOut(i) = 0;
///
//  Right edge point --- use boundary condition
//
        i = M;
	uOut(i) = 0;
//
//  Compute the vector z
//
//  Interior points: apply the operator
//

	i = 1;
	z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i));

	for(i = 2; i <= M-2; i++)
	{

		z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(i-1));
	}

	i = M-1;
	z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*( - 2.0*uIn(i) + uIn(i-1));

//
//  Add boundary conditions g0(t) and g1(t) to z
//  Both zeros for this problem
//

//
//  Given z, solve Ax = z
//

//  set up a_diag, b, c for the matrix A
	for (i=1; i <= M-1; i++)
	{
		a_diag(i-1) = 1.0 - (k/2.0)*(-a)/pow(h,2)*(-2.0);
		b(i-1) = - (k/2.0)*(-a)/pow(h,2)*(1.0);
		c(i-1) = - (k/2.0)*(-a)/pow(h,2)*(1.0);
	}

	TridiagSolver(a_diag, b, c, z, x);

//	set the output of the routine to the solution returned by the solver
	for(i = 1; i <= M-1; i++)
	{
		
		uOut(i) = x(i-1); 
	}
}

double maxNorm(DoubleArray1D& uIn, double h)
{
	long M = uIn.getSize();  // the whole vector
	long i;

	double result = 0.0;

	for (i = 0; i <= M; i++)
	{
		if (fabs(uIn(i)) >= result)
			result = fabs(uIn(i));
	}
	
	return result;
}

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
                e(i) = aIn(i) - d(i)*cIn(i-1);
        }

// forward substitution
        i = 0;
        y(i) = zIn(i);

        for(i = 1; i <= N-1; i++)
        {
                y(i) = zIn(i) - d(i)*y(i-1);
        }

// Backward substitution
        i = N-1;
        x(i) = y(i)/e(i);

        for(i = N-2; i >= 0; i--)
        {

                x(i) = (y(i) - cIn(i)*x(i+1))/e(i);
        }

//
        for(i = 0; i <= N-1; i++)
        {

                xOut(i) = x(i);
        }
}

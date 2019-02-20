//
//#####################################################################
//                            Assign6.cpp 
//#####################################################################
//
// This program computes the exact solution to 
// 
// u_t + a u_x = 0 
//
// for xA <= x <= xB, 0 <= t <= tFinal with 
// u(0,x) = u0(x) and u(t,x) being periodic in x 
// with period xB-xA (e.g. "periodic" boundary
// conditions). 
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
// Version: Feb. 10, 2010
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
double discreteL2Norm(DoubleArray1D&, double);

int main(int argc, char* argv[]) 
{
    double tFinal;              // final time

    double xA       = 0.0;      // left endpoint
    double xB       = 1.0;      // right endpoint
    double a        = -0.05;      // convection speed

    double lambda = 3.0/4.0;

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

    DoubleArray1D uExact;                       // Array to hold Exact solution values
    DoubleArray1D uN;                           // Array to hold computed solution values
    DoubleArray1D uNp1;                         // Array to hold computed next step solution values
    DoubleArray1D uErr;

    SinExactPeriodicSolution exactSoln(xA,xB,a);       // Exact solution class instance
    const char*  outputExactFileNamePrefix = "uExact";   // Exact output filename prefix.
    const char*  outputFileNamePrefix = "uComp";   // output filename prefix.
    int  outputFormat;                              // to specify output format


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

    dt = lambda * dx;
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
    uN.initialize(M+1);  // M+1 data points in our domain with M panels
    uNp1.initialize(M+1);  // M+1 data points in our domain with M panels
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
        uErr_norm = discreteL2Norm(uErr, dx);

        printf("Output %3ld Time = %10.5f  Timesteps taken : %ld \n", outputIndex, time, kStep);

        outputFormat = ProfileOutput::GNUPLOT;
        profileOut.output(uN, outputIndex, outputFileNamePrefix, outputFormat);
        profileOut.output(uExact, outputIndex, outputExactFileNamePrefix, outputFormat);
        printf("The discrete L2 norm of the error vector of length %3ld at time %5.2f is %10.5f \n", M+1, time, uErr_norm);
 
    }

    }


    return 0;
 }

//
// This routine advances the solution to u_t + a u_x = 0 one timestep using
// a method based upon a forward difference time derivative approximation and
// a one-sided difference spatial derivative approximation.
// 
// The solution u is assumed to be periodic. 
//
void advance(DoubleArray1D& uIn, DoubleArray1D& uOut, double h, double k, double a)
{
	long M = uIn.getSize() -1; // Number of panels = number of points - 1
	long i;
//
//  Left edge point --- use periodicity
//
	i = 0;
	uOut(i) = uIn(i) - (a*k)*(uIn(i+1)-uIn(M-1))/(2.0*h) + pow(a, 2)*pow(k,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(M-1))/(2.0*pow(h,2));
///
//  Right edge point --- use periodicity
//
        i = M;
	uOut(i) = uIn(i) - (a*k)*(uIn(1)-uIn(i-1))/(2.0*h) + pow(a, 2)*pow(k,2)*(uIn(1) - 2.0*uIn(i) + uIn(i-1))/(2.0*pow(h,2));
//
//  Interior points
//
	for(i = 1; i <= M-1; i++)
	{
		
		uOut(i) = uIn(i) - (a*k)*(uIn(i+1)-uIn(i-1))/(2.0*h) + pow(a, 2)*pow(k,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(i-1))/(2.0*pow(h,2));
	}
}

double discreteL2Norm(DoubleArray1D& uIn, double h)
{
	long M = uIn.getSize() -1;  // Number of panels = number of points - 1
	long i;

	double result = 0.0;

	for (i = 0; i <= M; i++)
	{
		result += pow(uIn(i), 2);	
	}
	
	result = sqrt(result * h);
	return result;
}

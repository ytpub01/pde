//
//#####################################################################
//                            Assign1.cpp 
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
//         uExactXXX.dat 
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
// Version: Jan. 31, 2010
// Author : Chris Anderson 
//#####################################################################
#include <iostream>
using namespace std;

#include <math.h>
#include <stdio.h>

#include "DoubleArray1D.h"         // 1D array class 
#include "ExactPeriodicSolution.h" // Class for the exact solution 

#include "ProfileOutput.h"         // Utility class for outputting data files in 
                                   // an appropriate format. 

int main(int argc, char* argv[]) 
{
    double tFinal;              // final time

    double xA       = 0.0;      // left endpoint
    double xB       = 1.0;      // right endpoint
    double a        = 1.0;      // convection speed

    long M;                     // number of panels
    double dx;                  // mesh spacing 

    double    time;             // total simulation time
    double      dt;             // timestep
    double  dtTemp;             // timestep for sub-stepping 

    long   outputCount;         // number of timsteps output
    double outputTimeIncrement; // time interval between output times
    double outputTime;          // next output time
    long   outputIndex;         // index of output steps

    DoubleArray1D uExact;                           // Array to hold exact solution values
    ExactPeriodicSolution exactSoln(xA,xB,a);       // Exact solution class instance
    const char*  outputFileNamePrefix = "uExact";   // output filename prefix.
    int  outputFormat;                              // to specify output format


    cout << " Enter number of grid panels M : " << endl;
    cin  >> M;
 
    cout << " Enter timestep size dt : " << endl;
    cin >> dt;

    cout <<  "Enter Final Time : " << endl;
    cin  >>  tFinal;

    cout <<  "Number of Output Times : " << endl;
    cin  >>  outputCount;

    //
    // Initialize spatial variables 
    //

    dx = (xB-xA)/double(M);

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
    exactSoln.evaluateExactSolution(uExact,0.0);

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
    profileOut.output(uExact, outputIndex, outputFileNamePrefix, outputFormat);
 
    printf("Output %3ld Time = %10.5f  Timesteps taken : %d \n", outputIndex, time,0);
    
    // Main time-stepping loop 

    int outputFlag       = 0;
    long kStep           = 0;
    int exitFlag         = 0;                   // XXXXX Correction 1/31/10
    while((time < tFinal)&&(exitFlag == 0))     // XXXXX Correction 1/31/10
    {

    // Simulation evolution 

    if(time + dt < outputTime)              // next step requires no output   
    {
        //
        // Advance numerical solution 
        //

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
        

        printf("Output %3ld Time = %10.5f  Timesteps taken : %ld \n", outputIndex, time, kStep);

        outputFormat = ProfileOutput::MATLAB;
        profileOut.output(uExact, outputIndex, outputFileNamePrefix, outputFormat);
 
    }

    }

    return 0;
 }


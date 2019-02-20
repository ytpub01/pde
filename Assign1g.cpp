//
//#####################################################################
//  Assign1g.cpp : This program uses the CAMgraphics classes and requires 
//  specific project settings and support files. See Assignment 1 
//  for details. 
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
// xA                                 xB
//   
// |---x---x---x---x---x---x---x---x---x
// 0   1   2   3   4   5               M
//
// 
//#####################################################################
// Created for Math 269B. 
// Version: Jan. 08, 2007
// Author: Chris Anderson 
//#####################################################################

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#include <math.h>
#include <stdio.h>

#include "DoubleArray1D.h"         // 1D array class 

#include "ExactPeriodicSolution.h" // Class for the exact solution 

#include "gprocess.h"             // CAMgraphics classes
#include "CAMmsvcDriver.h"        // Visual C++ Window Driver 


int CAMmain(int argc, char* argv[]) 
{
    long i;

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

    DoubleArray1D uExact;       //  Array to hold exact solution values
    DoubleArray1D xPlot;        //  Array to hold x values of grid for plotting
    DoubleArray1D yPlot;        //  Array to hold function values for plotting 

    ostringstream  titleStringStream; // String stream for plot titles
    char titleString[256];            // String for title

    ExactPeriodicSolution exactSoln(xA,xB,a); // Exact solution class instance

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
    // Graphics and plot variables initializion 
    //

    CAMgraphicsProcess G;           // declare a graphics process
    CAMmsvcDriver  Mdriver;         // declare an MSVC driver
    G.attachDriver(Mdriver);        // attach driver to process

    xPlot.initialize(M+1);
    yPlot.initialize(M+1);
    for(i = 0; i < M+1; i++)
    {
    xPlot(i) = xA + double(i)*dx;
    }

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
    // Plot (or output) the initial solution
    //

    double xMin = xA;
    double xMax = xB;
    double yMin = -1.0;
    double yMax =  3.0; 

    // capture initial solution 

    yPlot  = uExact;

    // create title 

    titleStringStream.str("");
    titleStringStream << "Time T = " << 0.0;
    strcpy(titleString,(titleStringStream.str()).c_str());

    // create plot 

    G.setAxisRange(xMin,xMax,yMin,yMax);
    G.labelsOn();
    G.title(titleString);
    G.plot(xPlot.getDataPointer(),yPlot.getDataPointer(),M+1);
    G.frame();

    printf("Output %3d Time = %10.5f  Timesteps taken : %d \n", outputIndex, time,0);

 
    // Main time-stepping loop 

    int outputFlag       = 0;
    long kStep           = 0;
    while(time < tFinal) 
    {

    // Simulation evolution 

    if(time + dt < outputTime)              // next step requires no output   
    {
        //
        // Advance numerical solution 
        //

        // update time and exact solution 
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
    }

    // Simulation output 

    if(outputFlag  == 1)
    {
        outputFlag  = 0;                    // reset flags and time markers
        outputTime += outputTimeIncrement;
        outputIndex++;      
        

        printf("Output %3d Time = %10.5f  Timesteps taken : %d \n", outputIndex, time, kStep);

        // capture solution 

        yPlot  = uExact;

        // create title for plot 

        titleStringStream.str("");
        titleStringStream << "Time T = " << time;
        strcpy(titleString,(titleStringStream.str()).c_str());

        // plot solution 

        G.setAxisRange(xMin,xMax,yMin,yMax);
        G.title(titleString);
        G.plot(xPlot.getDataPointer(),yPlot.getDataPointer(),M+1);
        G.labelsOn();
        G.frame();
    }

    }

    return 0;
 }

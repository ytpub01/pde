//
//#####################################################################
//                            Assign8.cpp 
//#####################################################################
//
// This program computes the exact solution to 
// 
// u_t + a (u_xx + u_yy) = 0 
//
// for xA <= x <= xB, yA <= y <= yB, and
// 0 <= t <= tFinal with 
// u(0,x,y) = u0(x,y), u(0,0,t) = g0(t), u(1,0,t) = g1(t)
// u(0,1,t) = g2(t)  u(1,1,t) = g3(t)   
// boundary conditions. 
//
// The spatial mesh is uniform with mPanel panels between
// x = xA and x = xB and nPanels between y = yA and y = yB.
// The array indices run from 0 to M as indicated below: 
// 
//
//               mPanel panels 
// xA                                  xB
//   
// |---x---x---x---x---x---x---x---x---x
// 0   1   2   3   4   5               mPanel  
//
//
//               nPanel panels 
// yA                                  yB
//   
// |---x---x---x---x---x---x---x---x---x
// 0   1   2   3   4   5               nPanel  
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
// Version: Mar. 3, 2010
// Author : Chris Anderson 
// Modified by : Yasser Taima
//#####################################################################
#include <iostream>
using namespace std;

#include <math.h>
#include <stdio.h>

#include "DoubleArray2D.h"         // 2D array class 
#include "DoubleArray1D.h"         // 1D array class 
#include "SinExactPeriodicSolution.h" // Class for the exact solution 

#include "ProfileOutput.h"         // Utility class for outputting data files in 
                                   // an appropriate format. 
#include "LaplaceOp2D_IB.h"              // 2D discrete Laplace operator

#include "outputToGNUplot.h"             // GNUplot output utility
#include "outputToMatlab.h"              // Matlab output utility

void advance2D(DoubleArray2D&, DoubleArray2D&, double, double, double, double, double, double, double, double);
double discreteL2Norm2D(const DoubleArray2D&, double, double);
double discreteL2Norm(DoubleArray1D&, double);
double maxNorm2D(DoubleArray2D&, double, double);
string composeFileName(long, const char*, int);
double sinPQ(double, double, double, double);

int main(int argc, char* argv[]) 
{
    double tFinal;              // final time

    double xA       = 0.0;      // left endpoint in the x direction
    double xB       = 1.0;      // right endpoint in the x direction
    double yA       = 0.0;      // left endpoint in the y direction
    double yB       = 1.0;      // right endpoint in the y direction
    double a        = -0.02;      // diffusion constant

    double p = 1.0;
    double q = 2.0;
    double pi = 4.0*atan(1.0);

    double myu;
    long mPanel;                // number of panels in the x direction
    long nPanel;                // number of panels in the y direction
    double hx;                  // mesh spacing in the x direction 
    double hy;                  // mesh spacing in the y direction 

    long i, j;
    double x, y;
    string s1;

    double    time;             // total simulation time
    double      dt;             // timestep
    double  dtTemp;             // timestep for sub-stepping 
    double  vErr_norm;

    long   outputCount;         // number of timsteps output
    double outputTimeIncrement; // time interval between output times
    double outputTime;          // next output time
    long   outputIndex;         // index of output steps

    DoubleArray2D vExact;       // Array to hold Exact solution values
    DoubleArray2D vN;           // Array to hold computed solution values
    DoubleArray2D vNp1;         // Array to hold computed next step solution values
    DoubleArray2D vErr;

    const char*  outputExactFileNamePrefix = "uExact";   // Exact output filename prefix.
    const char*  outputFileNamePrefix = "uComp";         // output filename prefix.
    int  outputFormat;                                   // to specify output format

    cout << " Enter number of grid panels in the x direction mPanel : " << endl;
    cin  >> mPanel;
 
    cout << " Enter number of grid panels in the y direction mPanel : " << endl;
    cin  >> nPanel;

 
    cout << " Enter myu : " << endl;
    cin  >> myu;
 
    cout <<  "Enter Final Time : " << endl;
    cin  >>  tFinal;

    cout <<  "Number of Output Times : " << endl;
    cin  >>  outputCount;

    //
    // Initialize spatial variables 
    //

    hx = (xB-xA)/double(mPanel);
    hy = (yB-yA)/double(nPanel);

    dt = myu * pow(hx, 2);

    vErr_norm = 0.0;

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

    vExact.initialize(mPanel+1, nPanel+1);  // mPanel+1 x nPanel+1 data points 
    vN.initialize(mPanel+1, nPanel+1);      
    vNp1.initialize(mPanel+1, nPanel+1);   
    vErr.initialize(mPanel+1, nPanel+1);

    //exact solution to begin
    for(j = 0; j <= nPanel; j++)
    {
    for(i = 0; i <= mPanel; i++)
    {
    y = j*hy + yA;
    x = i*hx + xA;
    vN(i,j) = sinPQ(x,y,p,q);
    vExact(i,j) = sinPQ(x,y,p,q);
    }}

    //
    // Instantiate output 
    //

    s1 =  composeFileName(outputIndex, outputFileNamePrefix, outputFormat);
    outputToMatlab(vN, s1.c_str());
    s1 =  composeFileName(outputIndex, outputExactFileNamePrefix, outputFormat);
    outputToMatlab(vExact, s1.c_str());

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
	advance2D(vN, vNp1, hx, hy, dt, a, xA, xB, yA, yB);
	vN = vNp1;

        // update time, stepcount, and exact solution 
        time += dt;
        kStep++;
	//advance exact solution
	for(j = 0; j <= nPanel; j++)
	{
	for(i = 0; i <= mPanel; i++)
	{
	y = j*hy + yA;
	x = i*hx + xA;
	vExact(i,j) = exp(-pow(pi,2)*(pow(p,2)+pow(q,2))*time*(-a))*sinPQ(x,y,p,q);
	}}
	//vN = vExact;
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

	advance2D(vN, vNp1, hx, hy, dt, a, xA, xB, yA, yB);
	vN = vNp1;

	// update time and exact solution 

        time += dtTemp;
        kStep++;
        //advance exact solution
	for(j = 0; j <= nPanel; j++)
        {
        for(i = 0; i <= mPanel; i++)
        {
        y = j*hy + yA;
        x = i*hx + xA;
	vExact(i,j) = exp(-pow(pi,2)*(pow(p,2)+pow(q,2))*time*(-a))*sinPQ(x,y,p,q);
        }}
        }
	if(fabs(time-tFinal) < 1.0e-10){exitFlag = 1;} // XXXXX Correction 1/31/10
	//vN = vExact;
    }

    // Simulation output 

    if(outputFlag  == 1)
    {
        outputFlag  = 0;                    // reset flags and time markers
        outputTime += outputTimeIncrement;
        outputIndex++;      
        
	vErr = vN - vExact;

        vErr_norm = discreteL2Norm2D(vErr, hx, hy);

        printf("Output %3ld Time = %10.5f  Timesteps taken : %ld \n", outputIndex, time, kStep);

	s1 =  composeFileName(outputIndex, outputFileNamePrefix, outputFormat);
	outputToMatlab(vN, s1.c_str());
	s1 =  composeFileName(outputIndex, outputExactFileNamePrefix, outputFormat);
	outputToMatlab(vExact, s1.c_str());

        printf("The discrete L2 norm of the error vector of length %3ld at time %5.2f is %10.5f \n", mPanel, time, vErr_norm);
 
    }

    }

    return 0;
 }

//
// This routine advances the solution to u_t + a u_xx = 0 one timestep using
// a method based upon a forward euler time derivative approximation and
// a two-sided difference spatial derivative approximation.
// 
//
void advance2D(DoubleArray2D& vIn, DoubleArray2D& vOut, double hx, double hy, double k, double a, double xA, double xB, double yA, double yB)
{
	long mPanel = vIn.getIndex1Size() -1; // Number of panels = number of points - 1
	long nPanel = vIn.getIndex2Size() -1; // Number of panels = number of points - 1
	long i, j;
	DoubleArray2D vInterm(mPanel+1, nPanel+1);
//
//  Southern edge --- use boundary condition
//
	j = 0;
	for (i = 0; i <= mPanel; i++)
		vIn(i,j) = 0.0;
//
//  Northern edge --- use boundary condition
//
        j = nPanel;
	for (i = 0; i <= mPanel; i++)
		vIn(i,j) = 0.0;
//
//  Eastern edge --- use boundary condition
//
	i = 0;
	for (j = 0; j <= nPanel; j++)
		vIn(i,j) = 0.0;
//
//  Western edge --- use boundary condition
//
	i = mPanel;
	for (j = 0; j <= nPanel; j++)
		vIn(i,j) = 0.0;
//
//	set the output of the routine

	LaplaceOp2D_IB laplaceOp(xA,xB,yA,yB);
	laplaceOp.apply(vIn, vInterm);

	vOut = vIn + (-a)*k*vInterm; 

}

double maxNorm2D(DoubleArray2D& vIn, double hx, double hy)
{
        long mPanel = vIn.getIndex1Size() -1;  // Number of panels = number of points - 1
        long nPanel = vIn.getIndex2Size() -1;  // Number of panels = number of points - 1
        long i, j;

        double result = 0.0;

        for (j = 0; j <= nPanel; j++)
	{
        	for (i = 0; i <= mPanel; i++)
        	{
                	if (fabs(vIn(i,j)) >= result)
                        	result = fabs(vIn(i,j));
        	}
	}

        return result;

}

double sinPQ(double x, double y, double p, double q)
{
    double pi =3.1415926535897932385;
    return sin(p*x*pi)*sin(q*y*pi);
}

double discreteL2Norm2D(const DoubleArray2D & V, double hx, double hy)
{
    long mPanel = V.getIndex1Size() + 1;
    long nPanel = V.getIndex2Size() + 1;

    double sum = 0.0;

    long i;
    long j;

    for(j = 0; j <= nPanel-2; j++)
    {
    for(i = 0; i <= mPanel-2; i++)
    {
    sum += V(i,j)*V(i,j)*hx*hy;
    }}
    return sqrt(sum);
}

string composeFileName(long outputIndex, const char* fileNamePrefix, int outputFormat)
{
    stringstream ss("");
    ss << fileNamePrefix;
    if     (outputIndex <= 9)  {ss << "00" << outputIndex << ".dat";}
    else if(outputIndex <= 99) {ss << "0" << outputIndex  << ".dat";}
    else                       {ss <<        outputIndex  << ".dat";}

    return ss.str();
}

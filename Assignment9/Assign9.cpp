//
//#####################################################################
//                            Assign9.cpp 
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
// Version: Mar. 10, 2010
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
void TridiagSolver(DoubleArray1D&, DoubleArray1D&, DoubleArray1D&, DoubleArray1D&, DoubleArray1D&);
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

    double lambda;
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

 
    cout << " Enter lambda or k (k <= 0.5): " << endl;
    cin  >> lambda;
 
    cout <<  "Enter Final Time: " << endl;
    cin  >>  tFinal;

    cout <<  "Number of Output Times: " << endl;
    cin  >>  outputCount;

    //
    // Initialize spatial variables 
    //

    hx = (xB-xA)/double(mPanel);
    hy = (yB-yA)/double(nPanel);

   if (lambda <= 0.5)
	dt = lambda;
   else
	dt = lambda * hx;
    
    printf("Choosing dt = %10.5f \n", dt);
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
	long i, p, q;	// p, q spatial indices
//  Arrays for the tridiagonal solver
        DoubleArray1D a_diag;                   // Array to hold diagonal values
        DoubleArray1D b;                        // Array to hold lower diagonal values
        DoubleArray1D c;                        // Array to hold upper diagonal values
        DoubleArray1D z;                        // Array to hold input vector
        DoubleArray1D x;                        // Array to hold solution

	DoubleArray2D vInterm(mPanel+1, nPanel+1); //intermediate structure to hold z in 2D
	DoubleArray2D vStar(mPanel+1, nPanel+1);

//  Southern edge --- use boundary condition
//
	q = 0;
	for (p = 0; p <= mPanel; p++)
	{
		vIn(p,q) = 0.0;
		vStar(p,q) = 0.0;
		vInterm(p,q) = 0.0;
	}
//
//  Northern edge --- use boundary condition
//
        q = nPanel;
	for (p = 0; p <= mPanel; p++)
	{
		vIn(p,q) = 0.0;
		vStar(p,q) = 0.0;
		vInterm(p,q) = 0.0;
	}
//
//  Eastern edge --- use boundary condition
//
	p = 0;
	for (q = 0; q <= nPanel; q++)
	{
		vIn(p,q) = 0.0;
		vStar(p,q) = 0.0;
		vInterm(p,q) = 0.0;
	}
//
//  Western edge --- use boundary condition
//

	p = mPanel;
	for (q = 0; q <= nPanel; q++)
	{
		vIn(p,q) = 0.0;
		vStar(p,q) = 0.0;
		vInterm(p,q) = 0.0;
	}
//
//	set the output of the routine

	// loop over p to get z in the q direction
	for(p=1; p<= mPanel-1; p++)	
	{
		// solve for q
		// set the z 1-D vector, a 1-D vector
		//
		//  Compute the vector z
		//
		//  Interior points: apply the operator
		//

	        i = 1;
	      	vInterm(p,i) = vIn(p,i) + (k/2.0)*(-a)/pow(hy,2)*(vIn(p,i+1) - 2.0*vIn(p,i)); 

	        for(i = 2; i <= nPanel-2; i++)
	                vInterm(p,i) = vIn(p,i) + (k/2.0)*(-a)/pow(hy,2)*(vIn(p,i+1) - 2.0*vIn(p,i) + vIn(p,i-1));

	        i = nPanel-1;
	        vInterm(p,i) = vIn(p,i) + (k/2.0)*(-a)/pow(hy,2)*( - 2.0*vIn(p,i) + vIn(p,i-1));
//
//  Add boundary conditions g0(t) and g1(t) to z
//  Both zeros for this problem
//
	}
//
//  Given z, solve Ax = z in the p direction, so loop over q
//

        a_diag.initialize(mPanel-1);  // N data points in the matrix
        b.initialize(mPanel-1);  // N-1 data points on the lower diagonal; b(0) is a dummy
        c.initialize(mPanel-1);  // N-1 data points on the lower diagonal; c(N-1) is a dummy
        z.initialize(mPanel-1);
        x.initialize(mPanel-1);

	for (q=1; q<= nPanel-1; q++)
	{
		for (i=1; i <= mPanel-1; i++)
			z(i-1) = vInterm(q,i); 
		//  set up a_diag, b, c for the matrix A1
        	for (i=1; i <= mPanel-1; i++)
        	{
	       		a_diag(i-1) = 1.0 - (k/2.0)*(-a)/pow(hx,2)*(-2.0);
        		b(i-1) = - (k/2.0)*(-a)/pow(hx,2)*(1.0);
        		c(i-1) = - (k/2.0)*(-a)/pow(hx,2)*(1.0); 
	        }

	        TridiagSolver(a_diag, b, c, z, x);

		// set v*(.,q) to the solution returned by the solver
        	for(i = 1; i <= mPanel-1; i++)
                	vStar(i,q) = x(i-1);
	}

	// loop over q to get z in the p direction
	for (q=1; q<= nPanel-1; q++)	
	{
		// solve for p
		// set the z 1-D vector, a 1-D vector
		//
		//  Compute the vector z
		//
		//  Interior points: apply the operator
		//

	        i = 1;
	      	vInterm(i,q) = vStar(i,q) + (k/2.0)*(-a)/pow(hx,2)*(vStar(i+1,q) - 2.0*vStar(i,q)); 

	        for(i = 2; i <= nPanel-2; i++)
	                vInterm(i,q) = vStar(i,q) + (k/2.0)*(-a)/pow(hx,2)*(vStar(i+1,q) - 2.0*vStar(i,q) + vStar(i-1,q));

	        i = nPanel-1;
	        vInterm(i,q) = vStar(i,q) + (k/2.0)*(-a)/pow(hx,2)*( - 2.0*vStar(i,q) + vStar(i-1,q));

//
//  Add boundary conditions g0(t) and g1(t) to z
//  Both zeros for this problem
//
	}
//
//  Given z, solve Ax = z in the q direction, so loop over p
//
        a_diag.initialize(nPanel-1);  // N data points in the matrix
        b.initialize(nPanel-1);  // N-1 data points on the lower diagonal; b(0) is a dummy
        c.initialize(nPanel-1);  // N-1 data points on the lower diagonal; c(N-1) is a dummy
        z.initialize(nPanel-1);
        x.initialize(nPanel-1);

	for (p=1; p<= mPanel-1; p++)
	{
		// set up z
		for (i=1; i <= nPanel-1; i++)
			z(i-1) = vInterm(p,i);
		//  set up a_diag, b, c for the matrix A2
	        for (i=1; i <= nPanel-1; i++)
        	{
		        a_diag(i-1) = 1.0 - (k/2.0)*(-a)/pow(hy,2)*(-2.0);
        		b(i-1) = - (k/2.0)*(-a)/pow(hy,2)*(1.0);
        		c(i-1) = - (k/2.0)*(-a)/pow(hy,2)*(1.0); 
        	}

        	TridiagSolver(a_diag, b, c, z, x);

		// set v*(p,.) to the solution returned by the solver
       		for(i = 1; i <= nPanel-1; i++)
               		vOut(p,i) = x(i-1);
	}	
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


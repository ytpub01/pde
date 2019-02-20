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
// Version: Mar. 2, 2010
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

void advance(DoubleArray1D&, DoubleArray1D&, double, double, double);
void advance2D(DoubleArray2D&, DoubleArray2D&, double, double, double, double, LaplaceOp2D_IB);
double maxNorm(DoubleArray1D&, double);
double discreteL2Norm(DoubleArray1D&, double);
double maxNorm2D(DoubleArray2D&, double, double);
void TridiagSolver(DoubleArray1D&, DoubleArray1D&, DoubleArray1D&, DoubleArray1D&, DoubleArray1D&);

int main(int argc, char* argv[]) 
{
    double pi =3.1415926535897932385;

    double tFinal;              // final time

    double xA       = 0.0;      // left endpoint in the x direction
    double xB       = 1.0;      // right endpoint in the x direction
    double yA       = 0.0;      // left endpoint in the y direction
    double yB       = 1.0;      // right endpoint in the y direction
    double a        = -0.02;      // diffusion constant

    double p = 1.0;
    double q = 2.0;

    LaplaceOp2D_IB laplaceOp(xA,xB,yA,yB);

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
    double  uErr_norm;

    long   outputCount;         // number of timsteps output
    double outputTimeIncrement; // time interval between output times
    double outputTime;          // next output time
    long   outputIndex;         // index of output steps

    DoubleArray1D uExact;       // Array to hold Exact solution values
    DoubleArray1D uN;           // Array to hold computed solution values
    DoubleArray1D uNp1;         // Array to hold computed next step solution values
    DoubleArray1D uErr;

    DoubleArray2D vExact;       // Array to hold Exact solution values
    DoubleArray2D vN;           // Array to hold computed solution values
    DoubleArray2D vNp1;         // Array to hold computed next step solution values
    DoubleArray2D vErr;

    SinExactPeriodicSolution exactSoln(xA,xB,a);         // Exact solution class instance
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
//    dt = hx;

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

    uExact.initialize(mPanel+1);  // mPanel+1 data points in our domain with mPanel panels
    uN.initialize(mPanel+1);      // mPanel+1 data points in our domain with mPanel panels
    uNp1.initialize(mPanel+1);    // mPanel+1 data points in our domain with mPanel panels
    uErr.initialize(mPanel+1);

    vExact.initialize(mPanel+1, nPanel+1);  // mPanel+1 data points in our domain with mPanel panels
    vN.initialize(mPanel+1, nPanel+1);      // mPanel+1 data points in our domain with mPanel panels
    vNp1.initialize(mPanel+1, nPanel+1);    // mPanel+1 data points in our domain with mPanel panels
    vErr.initialize(mPanel+1, nPanel+1);

    exactSoln.evaluateExactSolution(uExact,0.0);
    exactSoln.evaluateExactSolution(uN,0.0);
/*    for (j = 0; j <= nPanel; j++)
    {
    for (i = 0; i <= mPanel; i++)
    {
    vExact(i,j) = uExact(i);
    }}
    for (j = 0; j <= nPanel; j++)
    {
    for (i = 0; i <= mPanel; i++)
    {
    vN(i,j) = uN(i);
    }}
*/
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
//    outputFormat = ProfileOutput::GNUPLOT;
//    profileOut.output(uN, outputIndex, outputFileNamePrefix, outputFormat);
//    profileOut.output(uExact, outputIndex, outputExactFileNamePrefix, outputFormat);

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
        //
        // Advance numerical solution 
        //
/*	advance(uN, uNp1, hx, dt, a);
	uN = uNp1;
	for (j = 0; j <= nPanel; j++)
	{
	for (i = 0; i <= mPanel; i++)
	{
	vN(i,j) = uN(i);
	}}
*/
	advance2D(vN, vNp1, hx, hy, dt, a, laplaceOp);
	vN = vNp1;

        // update time, stepcount, and exact solution 
        time += dt;
        kStep++;
        exactSoln.evaluateExactSolution(uExact,time); 
	for (j = 0; j <= nPanel; j++)
	{
	for (i = 0; i <= mPanel; i++)
	{
	y = j*hy + yA;
	x = i*hx + xA;
//	vExact(i,j) = uExact(i);
	vExact(i,j) = exp(-pow(pi,2)*(pow(p,2)+pow(q,2))*(-a)*time)*sinPQ(x,y,p,q);
	}}

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

/*	for (j = 0; j <= nPanel; j++)
	{
		for (i = 0; i <= mPanel; i++)
			{
				uN(i) = vN(i, j);
			}
        	advance(uN, uNp1, hx, dtTemp, a);
		uN = uNp1;
                for (i = 0; i <= mPanel; i++)
                        {
                                vN(i,j) = uN(i);
                        }
	}
*/
	advance2D(vN, vNp1, hx, hy, dt, a, laplaceOp);
	vN = vNp1;

	// update time and exact solution 

        time += dtTemp;
        kStep++;
        exactSoln.evaluateExactSolution(uExact,time);
	for (j = 0; j <= nPanel; j++)
	{
	for (i = 0; i <= mPanel; i++)
	{
	y = j*hy + yA;
	x = i*hx + xA;
//	vExact(i,j) = uExact(i);
	vExact(i,j) = exp(-pow(pi,2)*(pow(p,2)+pow(q,2))*(-a)*time)*sinPQ(x,y,p,q);
	}}

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
	vErr = vN - vExact;

//        uErr_norm = discreteL2Norm(uErr, hx);
        uErr_norm = discreteL2Norm2D(vErr, hx, hy);
//        uErr_norm = maxNorm(uErr, hx);
//        uErr_norm = maxNorm2D(vErr, hx, hy);

        printf("Output %3ld Time = %10.5f  Timesteps taken : %ld \n", outputIndex, time, kStep);

//        outputFormat = ProfileOutput::GNUPLOT;
//        profileOut.output(uN, outputIndex, outputFileNamePrefix, outputFormat);
//        profileOut.output(uExact, outputIndex, outputExactFileNamePrefix, outputFormat);

	s1 =  composeFileName(outputIndex, outputFileNamePrefix, outputFormat);
	outputToMatlab(vN, s1.c_str());
	s1 =  composeFileName(outputIndex, outputExactFileNamePrefix, outputFormat);
	outputToMatlab(vExact, s1.c_str());

        printf("The max norm of the error vector of length %3ld at time %5.2f is %10.5f \n", mPanel, time, uErr_norm);
 
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
void advance2D(DoubleArray2D& vIn, DoubleArray2D& vOut, double hx, double hy, double k, double a, LaplaceOp2D_IB laplaceOp)
{
	long mPanel = vIn.getIndex1Size() -1; // Number of panels = number of points - 1
	long nPanel = vIn.getIndex2Size() -1; // Number of panels = number of points - 1
	long i, j;

//  Arrays for the tridiagonal solver
	DoubleArray2D vInterm;
	DoubleArray2D a_diag;                   // Array to hold diagonal values
	DoubleArray2D b;                        // Array to hold lower diagonal values
	DoubleArray2D c;                        // Array to hold upper diagonal values
	DoubleArray2D z;                        // Array to hold input vector
	DoubleArray2D x;                        // Array to hold solution
	
	a_diag.initialize(mPanel-1, nPanel-1);  // N data points in the matrix
	b.initialize(mPanel-1, nPanel-1);  // N-1 data points on the lower diagonal; b(0) is a dummy
	c.initialize(mPanel-1, nPanel-1);  // N-1 data points on the lower diagonal; c(N-1) is a dummy
	z.initialize(mPanel-1, nPanel-1); 
	x.initialize(mPanel-1, nPanel-1); 
	vInterm.initialize(mPanel+1, nPanel+1);
//
//  Southern edge --- use boundary condition
//
	j = 0;
	for (i = 0; i <= mPanel; i++)
		vOut(i,j) = 0.0;
//
//  Northern edge --- use boundary condition
//
        j = nPanel;
	for (i = 0; i <= mPanel; i++)
		vOut(i,j) = 0.0;
//
//  Eastern edge --- use boundary condition
//
	i = 0;
	for (j = 0; j <= nPanel; j++)
		vOut(i,j) = 0.0;
//
//  Western edge --- use boundary condition
//
	i = mPanel;
	for (j = 0; j <= nPanel; j++)
		vOut(i,j) = 0.0;
//
//
//  Compute the vector z
//
//  Interior points: apply the operator
//

/*	i = 1;
	// F Euler
	z(i-1) = uIn(i) + k*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i));

	// B Euler
	//z(i-1) = uIn(i);

	// Trapezoidal
        //z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i));

	for(i = 2; i <= mPanel-2; i++)
	{
		// F Euler
		z(i-1) = uIn(i) + k*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(i-1));

		// B Euler
		//z(i-1) = uIn(i);

		// Trapezoidal
                //z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(i-1));
	}

	i = mPanel-1;
	// F Euler
	z(i-1) = uIn(i) + k*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(i-1));
	
	// B Euler
	//z(i-1) = uIn(i);
	
	// Trapezoidal
        //z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*( - 2.0*uIn(i) + uIn(i-1));

//
//  Add boundary conditions g0(t) and g1(t) to z
//  Both zeros for this problem
//

//
//  Given z, solve Ax = z
//

//  set up a_diag, b, c for the matrix A
	for (i=1; i <= mPanel-1; i++)
	{
		// F Euler 
		a_diag(i-1) = 1.0;
		b(i-1) = 0.0;
		c(i-1) = 0.0;

		// B Euler 
		//a_diag(i-1) = 1.0 - k*(-a)/pow(h,2)*(-2.0);
		//b(i-1) = - k*(-a)/pow(h,2)*(1.0);
		//c(i-1) = - k*(-a)/pow(h,2)*(1.0); 

		// Trapezoidal
		//a_diag(i-1) = 1.0 - (k/2.0)*(-a)/pow(h,2)*(-2.0);
		//b(i-1) = - (k/2.0)*(-a)/pow(h,2)*(1.0);
		//c(i-1) = - (k/2.0)*(-a)/pow(h,2)*(1.0); 
	}

	TridiagSolver(a_diag, b, c, z, x);

//	set the output of the routine to the solution returned by the solver
	for(i = 1; i <= mPanel-1; i++)
	{
		
		uOut(i) = x(i-1); 
	}
*/
//	set the output of the routine
	laplaceOp.apply(vIn, vInterm);
	vOut = vIn + k*(-a)*vInterm;
}

void advance(DoubleArray1D& uIn, DoubleArray1D& uOut, double h, double k, double a)
{
	long mPanel = uIn.getSize() -1; // Number of panels = number of points - 1
	long i;

//  Arrays for the tridiagonal solver
	DoubleArray1D a_diag;                   // Array to hold diagonal values
	DoubleArray1D b;                        // Array to hold lower diagonal values
	DoubleArray1D c;                        // Array to hold upper diagonal values
	DoubleArray1D z;                        // Array to hold input vector
	DoubleArray1D x;                        // Array to hold solution
	
	a_diag.initialize(mPanel-1);  // N data points in the matrix
	b.initialize(mPanel-1);  // N-1 data points on the lower diagonal; b(0) is a dummy
	c.initialize(mPanel-1);  // N-1 data points on the lower diagonal; c(N-1) is a dummy
	z.initialize(mPanel-1); 
	x.initialize(mPanel-1); 
//
//  Left edge point --- use boundary condition
//
	i = 0;
	uOut(i) = 0;
///
//  Right edge point --- use boundary condition
//
        i = mPanel;
	uOut(i) = 0;
//
//  Compute the vector z
//
//  Interior points: apply the operator
//

	i = 1;
	// F Euler
	z(i-1) = uIn(i) + k*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i));

	// B Euler
	//z(i-1) = uIn(i);

	// Trapezoidal
        //z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i));

	for(i = 2; i <= mPanel-2; i++)
	{
		// F Euler
		z(i-1) = uIn(i) + k*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(i-1));

		// B Euler
		//z(i-1) = uIn(i);

		// Trapezoidal
                //z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(i-1));
	}

	i = mPanel-1;
	// F Euler
	z(i-1) = uIn(i) + k*(-a)/pow(h,2)*(uIn(i+1) - 2.0*uIn(i) + uIn(i-1));
	
	// B Euler
	//z(i-1) = uIn(i);
	
	// Trapezoidal
        //z(i-1) = uIn(i) + (k/2.0)*(-a)/pow(h,2)*( - 2.0*uIn(i) + uIn(i-1));

//
//  Add boundary conditions g0(t) and g1(t) to z
//  Both zeros for this problem
//

//
//  Given z, solve Ax = z
//

//  set up a_diag, b, c for the matrix A
	for (i=1; i <= mPanel-1; i++)
	{
		// F Euler 
		a_diag(i-1) = 1.0;
		b(i-1) = 0.0;
		c(i-1) = 0.0;

		// B Euler 
		//a_diag(i-1) = 1.0 - k*(-a)/pow(h,2)*(-2.0);
		//b(i-1) = - k*(-a)/pow(h,2)*(1.0);
		//c(i-1) = - k*(-a)/pow(h,2)*(1.0); 

		// Trapezoidal
		//a_diag(i-1) = 1.0 - (k/2.0)*(-a)/pow(h,2)*(-2.0);
		//b(i-1) = - (k/2.0)*(-a)/pow(h,2)*(1.0);
		//c(i-1) = - (k/2.0)*(-a)/pow(h,2)*(1.0); 
	}

	TridiagSolver(a_diag, b, c, z, x);

//	set the output of the routine to the solution returned by the solver
	for(i = 1; i <= mPanel-1; i++)
	{
		
		uOut(i) = x(i-1); 
	}
}

double discreteL2Norm(DoubleArray1D& uIn, double h)
{
        long mPanel = uIn.getSize() -1;  // Number of panels = number of points - 1
        long i;

        double result = 0.0;

        for (i = 0; i <= mPanel; i++)
        {
                result += pow(uIn(i), 2);
        }

        result = sqrt(result * h);
        return result;
}

double maxNorm(DoubleArray1D& uIn, double h)
{
        long mPanel = uIn.getSize() -1;  // Number of panels = number of points - 1
        long i;

        double result = 0.0;

        for (i = 0; i <= mPanel; i++)
        {
                if (fabs(uIn(i)) >= result)
                        result = fabs(uIn(i));
        }

        return result;
}

double maxNorm2D(DoubleArray2D& vIn, double hx, double hy)
{
        long mPanel = vIn.getIndex1Size() -1;  // Number of panels = number of points - 1
        long nPanel = vIn.getIndex2Size() -1;  // Number of panels = number of points - 1
        long i, j;

        double result = 0.0;

        for (j = 0; j <= mPanel; j++)
	{
        	for (i = 0; i <= mPanel; i++)
        	{
                	if (fabs(vIn(i,j)) >= result)
                        	result = fabs(vIn(i,j));
        	}
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

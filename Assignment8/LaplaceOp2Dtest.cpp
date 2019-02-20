//
//#####################################################################
//                     LaplaceOp2DTest.cpp 
//
// =======    Only interior data values are stored   ======
//
//#####################################################################
//
// This program computes tests the 2D discrete Laplace operator
// for discrete functions defined on the domain [xA, xB]x[yA, yB]
// with a uniform mesh with mPanel panels in the x-direction and
// nPanel panels in the y-direction.  
//
// In this implementation, only the interior values are stored in a
// single DoubleArray2D. The boundary values are stored in separate
// arrays bdXA, bdXB, bdYA, bdYB.
// 
// For interior values:
// the i array index (x-direction) runs from 0 to mPanel-2.
// the j array index (y-dicretion) runs from 0 to nPanel-2.
// 
// This program also demonstrates the use of utility functions
// for creating files containing 2D data that can be read and
// plotted by GNUplot or Matlab.
//
// Currently the domain size is fixed to be [0, 1]x[0, 1]
//
//#####################################################################
// Initial Construction : Feb. 23, 2010
// 2/28/10: Corrected boundary indexing in Vplot packing
//#####################################################################


#include <iostream>
using namespace std;

#include <math.h>
#include <stdio.h>

#include "DoubleArray2D.h"                // 2D array class 
#include "DoubleArray1D.h"                // 1D array class 

#include "LaplaceOp2D.h"                  // 2D discrete Laplace operator

#include "outputToGNUplot.h"             // GNUplot output utility
#include "outputToMatlab.h"              // Matlab output utility


// Test function 

double sinPQ(double x, double y, double p, double q)
{
    double pi =3.1415926535897932385;
    return sin(p*x*pi)*sin(q*y*pi);
}
//
// This routine computes the discrete L2 Norm of the interior
// values of a function associated with a uniform rectangular grid with mesh
// size hx in the x-direction and hy in the y-direction.
//
// The input DoubleArray2D V is assumed to be of size mPanel-1 by nPanel-1.
//
// 
double discreteL2Norm(const DoubleArray2D & V, double hx, double hy)
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

int main(int argc, char* argv[]) 
{

    double xA       = 0.0;      // left endpoint
    double xB       = 1.0;      // right endpoint
    double yA       = 0.0;      // bottom endpoint
    double yB       = 1.0;      // top endpoint 
    double b        = 1.0;      // diffusion coefficient

    long mPanel     = 10;       // number of panels in x-direction
    double hx;                  // mesh spacing in x-direction

    long nPanel     = 20;       // number of panels in y-direction
    double hy;                  // mesh spacing in y-direction

    hx = (xB-xA)/double(mPanel);
    hy = (yB-yA)/double(nPanel);

//
//  Check the discrete Laplace operator by applying it to 
//  sin(p*pi*x)*sin(q*pi*y) -- the discrete function associated
//  with this function is an eigenvector of the discrete operator.
//
//  The only error in the computation should be that due to 
//  numerical imprecision.
//

    double p = 2.0;
    double q = 3.0;

    double pi = 4.0*atan(1.0);
    
    double lambdaPQ;  // eigenvalue for the (p,q)th eigenvector
    
    lambdaPQ  = -4.0*(sin(p*pi*hx*0.5)*sin(p*pi*hx*0.5))/(hx*hx)
                -4.0*(sin(q*pi*hy*0.5)*sin(q*pi*hy*0.5))/(hy*hy);


    DoubleArray2D Vin(mPanel-1,nPanel-1);
    DoubleArray2D Vout(mPanel-1,nPanel-1);

    DoubleArray1D bdXA(nPanel-1); // bdry data at x = xA
    DoubleArray1D bdXB(nPanel-1); // bdry data at x = xB

    DoubleArray1D bdYA(mPanel-1); // bdry data at y = yA
    DoubleArray1D bdYB(mPanel-1); // bdry data at y = yB

    //
    // Specify interior values
    //
    long i; long j;
    double x; 
    double y;

    for(j = 0; j <= nPanel-2; j++)
    {
    for(i = 0; i <= mPanel-2; i++)
    {
    x = (i+1)*hx + xA;
    y = (j+1)*hy + yA;
    Vin(i,j) = sinPQ(x,y,p,q);
    }}
    //
    // Specify homogeneous boundary values
    //
    bdXA.setToValue(0.0);
    bdXB.setToValue(0.0);
    bdYA.setToValue(0.0);
    bdYB.setToValue(0.0);

    LaplaceOp2D laplaceOp(xA,xB,yA,yB);
//
// Apply the Laplace operator using the specified boundary values
//
    laplaceOp.apply(Vin,bdXA, bdXB, bdYA, bdYB, Vout);
//
//  Check results
//
    DoubleArray2D Verr(mPanel-1,nPanel-1);
   
    Verr  = Vout - lambdaPQ*Vin;

    double VNorm    = discreteL2Norm(Vin,hx,hy);
    double VerrNorm = discreteL2Norm(Verr,hx,hy);

    cout << " Norm of input function :  " << VNorm << endl;
    cout << " Norm of Error in Discrete Laplace Operator: " << VerrNorm << endl;
//
//  Create array containing interior and boundary values for output
//
    DoubleArray2D Vplot(mPanel+1,nPanel+1);

    for(j = 0; j <= nPanel-2; j++)
    {
    for(i = 0; i <= mPanel-2; i++)
    {
    Vplot(i+1,j+1) = Vout(i,j);
    }}

    for(j = 0; j <= nPanel-2; j++)
    {
    i = 0;
    Vplot(i,j+1) = bdXA(j);
    i = mPanel;
    Vplot(i,j+1) = bdXB(j);
    }

    for(i = 0; i <= mPanel-2; i++)
    {
    j = 0;
    Vplot(i+1,j) = bdYA(i);
    j = nPanel;
    Vplot(i+1,j) = bdYB(i);
    }

    Vplot(0,0)           = 0.0;
    Vplot(0,nPanel)      = 0.0;
    Vplot(0,nPanel)      = 0.0;
    Vplot(mPanel,nPanel) = 0.0;
//
//  Save the input data to a file that can be read by GNUplot's splot 
//  command.
//
//  Set up coordinates locations for plot
//
    DoubleArray1D Xplot(mPanel+1);  // x coordinates
    DoubleArray1D Yplot(nPanel+1);  // y coordinates
   
    for(i = 0; i <= mPanel; i++)
    {
    Xplot(i) = i*hx + xA;
    }
    for(j = 0; j <= nPanel; j++)
    {
    Yplot(j) = j*hy + yA;
    }

//
//  Call the GNUplot output routine 
//
    outputToGNUplot(Xplot, Yplot, Vplot, "GNU_test2D.dat");

/* 
//    To view the plot use the commands

       >> set pm3d
       >> splot 'GNU_test2D.dat' w l

       or if your GNUplot doesn't support pm3d, just use
   
       >> splot 'GNU_test2D.dat' w l 
*/
//
//  Call the Matlab output routine
//
    outputToMatlab(Vplot, "matlab_test2D.dat");

/*
//  To view the plot, go to the directory containing the
    output files and use the commands
 
       >> load 'matlab_test2D.dat'
       >> surf(matlab_test2d)

    Note: if after loading, surf doesn't find matlab_test2d, then execute
       >> who 
    to discover how Matlab has named the file. 
*/

    return 0;
 }

 

#include <iostream>
using namespace std;

#include <math.h>          // Using C math functions

#include "DoubleArray1D.h" // Using DoubleArray1D class
//
// Routine 
//
//               applyBackwardDiffOp
//
// Apply a backward difference operator to the function
// values associated with a periodic grid with mesh
// spacing h and a number of panels M.  
//
// Input : 
//
// h     : The mesh size associated with a uniform periodic
//         grid. 
// 
// u     : DoubleArray1D instance of size M+1 whose elements 
//         correspond the the values of a periodic grid function. 
//         It is assumed that u(M) = u(0)
//
// Output :
//
// Returns a DoubleArray1D instance Du with values that are 
// obtained by applying a backward difference operator to u :
//
// (Du)(i)  = (u(i) - u(i-1)/h for i = 0..M (*)
// 
// where periodicity is used to determine u(-1), e.g. u(-1)= u(M-1) 
//
// The vector Du is periodic, e.g. Du(0) = Du(M)
//
//
DoubleArray1D applyBackwardDiffOp(double h, DoubleArray1D& u)
{
    long M = u.getSize()- 1;   // M = number of panels
    DoubleArray1D Du(M+1);     // return argument
    
    Du(0) = (u(0)-u(M-1))/h;   // Difference at i = 0, using periodicity
                               // to determine u(-1) value. 
    
    long i;                    // Backward difference at the remaining points. 
    for(i = 1; i <= M; i++)
    {
        Du(i) = (u(i) - u(i-1))/h;
    }
    return Du;
}

//
// The evaluateCos routine evaluates cos(x) at the positions 
// specified by input DoubleArray1D vector x
//
DoubleArray1D evaluateCos(DoubleArray1D& x)
{
    long M = x.getSize()-1;   // size of input array
    DoubleArray1D cosU(M+1);  // return argument

    long i;
    for(i = 0; i <= M; i++) 
    {
    cosU(i) = cos(x(i));
    } 

    return cosU;
}

//
// This evaluateSin routine evaluates sin(x) at the positions 
// specified by input DoubleArray1D vector x. 
//
DoubleArray1D evaluateSin(DoubleArray1D& x)
{
    long M = x.getSize()-1;   // size of input array
    DoubleArray1D sinU(M+1);  // return argument

    long i;
    for(i = 0; i <= M; i++) 
    {
    sinU(i) = sin(x(i));
    } 

    return sinU;
}

// 
// Code copied from Assignment 1 to evaluate the 
// discrete 2 norm of a DoubleArray1D instance. 
//
double periodic2Norm(DoubleArray1D& u, double h)
{
    long M = u.getSize()-1;
    double sum = 0.0;
    long i;
    for(i = 0; i <= M-1; i++)
    {
    sum += u(i)*u(i)*h;
    }
    sum = sqrt(sum);
    return sum;
}

//
// Routine maxNorm implements the L infinity (max norm) of the
// vector u
//
double maxNorm(DoubleArray1D& u)
{
    long M = u.getSize()-1; // get panel count
    
    double maxU = 0.0;
    for(long i = 0; i <= M; i++)
    {
    maxU = (fabs(u(i)) > maxU ) ? fabs(u(i)) : maxU;
    }
    return maxU;
}

//
// Main code
//
int main()
{
    double pi = 4.0*atan(1.0);
    
    double xMin = 0.0;        // Create an interval over which sin is periodic
    double xMax = 4.0*pi;
   
    long M;
    cout << " Enter number of panels: ";
    cin  >> M;
   
    double h = (xMax-xMin)/(double)M;
    long i;
   
    DoubleArray1D x(M+1);       // for x locations
    DoubleArray1D u(M+1);       // for test function
    DoubleArray1D uComp(M+1);   // for computed result
    DoubleArray1D uExact(M+1);  // for exact result
    DoubleArray1D uErr(M+1);    // for error 
    
    for(i = 0; i <= M; i++)
    {
        x(i) = xMin + h*(double)i;
    }
    
    u      = evaluateSin(x);
    uExact = evaluateCos(x);
    
    uComp  =  applyBackwardDiffOp(h, u);
    
    uErr = uComp - uExact;
    
    double uNorm      = periodic2Norm(uExact, h);
    double errNorm    = periodic2Norm(uErr, h);
    double maxErrNorm = maxNorm(uErr);
    
    //
    // Set the number of digits for output (this does not set
    // the precision of the calculation, only the number of digits
    // that are displayed). 
    //
    cout.precision(15);
    
    cout << " Absolute error L2 norm  " << errNorm << endl;
    cout << " Relative error L2 norm  " << errNorm/uNorm << endl;
    cout << " Absolute error Max norm " << maxErrNorm << endl;

    return 0;
}

 

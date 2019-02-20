#ifndef __SinExactPeriodicSolution__
#define __SinExactPeriodicSolution__

#include "DoubleArray1D.h"
#include <math.h>
//
// Member functions of the SinExactPeriodicSolution class returns
// the exact solution to u_t + a u_x = 0 
// for x in the interval [xA, xB]. 
//
// The initial data implemented is the function sin((coeff pi (x-xA))/(xB-xA)).
// 
class SinExactPeriodicSolution
{
public :

    //
    // Constructor : Specify the interval [xA, xB] and the
    //               convection speed a. 
    //

    SinExactPeriodicSolution(double xA, double xB, double a)
    {
        this->xA        = xA;
        this->xB        = xB;
        this->xInterval = xB-xA;
        this->a         = a;
    }

    //
    // Given x in the interval [xA,xB] this function returns
    // the exact solution at time t. 
    //
    double evaluateExactSolution(double x, double t)
    {
        double pStar;

        // determine relative location in the
        // first periodic interval 
 
        // evaluate initial data at this location
	
	pStar = exp(-pow(1.0,2)*pow((3.1415926535897932385),2)*(-a)*t)*initialFunction(x, 1.0)
		- exp(-pow(3.0,2)*pow((3.1415926535897932385),2)*(-a)*t)*initialFunction(x, 3.0); 
//	pStar = exp(-pow(2.0,2)*pow((3.1415926535897932385),2)*(-a)*t)*initialFunction(x, 2.0);
        
	return pStar;
    }

    //
    // This routine samples the exact solution at time t at
    // the nodes of an equispaced grid whose number of panels
    // is one less than the size of the input array u. 
    //
    void evaluateExactSolution(DoubleArray1D& u, double t)
    {
        long mPanel = u.getSize() - 1;
        double dx   = (xB - xA)/double(mPanel);
        double xp;
    
        long i;

        for(i = 0; i < mPanel+1; i++)
        {
        xp = xA + double(i)*dx;
        u(i) = evaluateExactSolution(xp,t);
        }
        return;
    }

    //
    // The initial function is sin((coeff* pi (x-xA))/(xB-xA)),
    // The specification is with respect to the periodic interval 
    // [0, xInterval].
    //

    double initialFunction(double pStar, double coeff)
    {
        double s = pStar/(xInterval);   // normalize to [0,1]
        return sin(coeff*3.1415926535897932385*s);
    }

    double xA;        // periodic domain bounds
    double xB;
    double xInterval; // periodic size
    double a;         // convective wave speed
};

#endif

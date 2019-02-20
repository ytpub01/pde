#ifndef __ExactPeriodicSolution__
#define __ExactPeriodicSolution__

#include "DoubleArray1D.h"
#include <math.h>
//
// Member functions of the ExactPeriodicSolution class returns
// the exact solution to u_t + a u_x = 0 
// for x in the interval [xA, xB]. 
//
// The initial data implemented is a parabolic "hump"
//  centered in the interval. 
// 
class ExactPeriodicSolution
{
public :

    //
    // Constructor : Specify the interval [xA, xB] and the
    //               convection speed a. 
    //

    ExactPeriodicSolution(double xA, double xB, double a)
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

        pStar = x - a*t;

        // determine relative location in the
        // first periodic interval 

        pStar = pStar - xA;
        pStar = pStar - xInterval*floor(pStar/xInterval);
 
        // evaluate initial data at this location 

        return initialFunction(pStar);
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
    // The initial function is a parabolic hump centered in the
    // the middle of the domain. The specification is with respect
    // to the periodic interval [0, xInterval].
    //

    double initialFunction(double pStar)
    {
        double s = pStar/(xInterval);   // normalize to [0,1]
        if(s < 0.25) return 0.0;
        if(s > 0.75) return 0.0;
        return  16.0*(s-.25)*(.75-s);
    }

    double xA;        // periodic domain bounds
    double xB;
    double xInterval; // periodic size
    double a;         // convective wave speed
};

#endif


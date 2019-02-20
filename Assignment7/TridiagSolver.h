#ifndef __TridiagSolver__
#define __TridiagSolver__

#include "DoubleArray1D.h"
#include <math.h>
//
// Member functions of the TridiagSolver class returns
// the exact solution to Ax = z 
// for A an invertible tridiagonal matrix, z a vector in R^N. 
//
class TridiagSolver
{
public :

    //
    // Constructor : Specify the three vectors for the diagonal
    // upper diagonal and lower diagonal elements of the matrix A
    // as well as the vector z
    //

    TridiagSolver(DoubleArray1D& aIn, DoubleArray1D& bIn, DoubleArray1D& cIn)
    {
        this->aIn        = aIn;
        this->bIn        = bIn;
        this->cIn        = cIn;
    }

    //
    // Given z this function returns
    // the exact solution of Ax = z. 
    //
    DoubleArray1D& evaluateExactSolution(DoubleArray1D& z)
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
//        printf("d = %10.5f  Timestep : %ld \n", d(i), i);
                e(i) = aIn(i) - d(i)*cIn(i-1);
//        printf("e = %10.5f  Timestep : %ld \n", e(i), i);
        }

// forward substitution
        i = 0;
        y(i) = zIn(i);

        for(i = 1; i <= N-1; i++)
        {
                y(i) = zIn(i) - d(i)*y(i-1);
//        printf("y = %10.5f  Timestep : %ld \n", y(i), i);
        }

// Backward substitution
        i = N-1;
        x(i) = y(i)/e(i);
//        printf("x = %10.5f  Timestep : %ld \n", x(i), i);

        for(i = N-2; i >= 0; i--)
        {

                x(i) = (y(i) - cIn(i)*x(i+1))/e(i);
//        printf("x = %10.5f  Timestep: %ld \n", x(i), i);
        }

//
	return x;

    }

#endif

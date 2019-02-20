//
//######################################################################
//
// GNUplotArray2Dtest.cpp : a test program that illustrates the construction
//                     of a gnuplot readable data file that contains the
//                     the values of a scalar function of two variables
//                     sampled on a cartesian product mesh, e.g. at the
//                     points (X(i),Y(j)) for i = 0 ... mPanels and 
//                     j = 0 .. nPanels;
//                
// 
// Math 270C                                                  03/05/04
// Math 269B                                                  02/27/07
//                                 Fixed X,Y initialization   03/06/07
//######################################################################
//
//
#include "outputToGNUplot.h"

#include "DoubleArray1D.h"
#include "DoubleArray2D.h"
#include <math.h>

#include <iostream>
using namespace std;

double f(double x, double y)
{
	return cos(x*y)*sin(x*x + y*y);
}

int main()
{
//
// Specify the output file. If the file contains the DOS directory 
// characther \ then you'll need to use \\ in the string 
// to designate it. e.g. "H:\\CRA\\270c.1.04w\\CPPcodes\\GNUplotTest2D.dat"
// 
// The default output location is the directory in which the executable is 
// located. 
//
	char* outputFileName = "GNUplotTest2D.dat";

	long mPanel = 50;
	long nPanel = 50;

	double xMin = -2.0;
	double xMax =  2.0;

	double yMin = -4.0;
	double yMax =  4.0;

	double dx   = (xMax - xMin)/double(mPanel);
	double dy   = (yMax - yMin)/double(nPanel);
//
//  Set up evaluation points
//
	DoubleArray1D X(mPanel+1);
	DoubleArray1D Y(nPanel+1);

	long i; long j;

	for(i = 0; i <= mPanel; i++)
	{
		X(i) = xMin + dx*(i);
	}

	for(j = 0; j <= nPanel; j++)
	{
		Y(j) = yMin + dy*(j);
	}

	DoubleArray2D F(mPanel+1,nPanel+1);
//
//  Sample the function to be plotted
//
	for(i = 0; i <= mPanel; i++)
	{
	for(j = 0; j <= nPanel; j++)
	{
		F(i,j)= f(X(i),Y(j));
	}
	}

	outputToGNUplot(X, Y, F,outputFileName);

	return 0;
}

 

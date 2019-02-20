//
//######################################################################
//
// outputToGNUplot.cpp : a routine that constructs of a gnuplot readable 
// data file that contains the the values of a scalar function of two 
// variables sampled on a cartesian product mesh, e.g. at the
// points (X(i),Y(j)) for i = 0 ... mPanels and j = 0 .. nPanels;
//                
// 
// Math 270C                                                   03/05/04
// Math 269B                                                   02/23/10
//######################################################################
//
//

#ifndef __OutputToGNUplot__
#define __OutputToGNUplot__

#include "DoubleArray1D.h"
#include "DoubleArray2D.h"

void outputToGNUplot(DoubleArray1D& Xcoord, DoubleArray1D&Ycoord, 
DoubleArray2D& F, const char* fileName);

#endif
 

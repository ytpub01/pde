//
//######################################################################
//
// outputToMatlab.h : header for a routine that constructs a Matlab readable 
// data file that contains the values of a scalar function of two 
// variables sampled on a cartesian product mesh, e.g. at the
// grid points (i,j) for i = 0 ... mPanels and j = 0 .. nPanels;
//                
// 
// Math 269B                                                   02/23/10
//######################################################################
//
//

#ifndef __OutputToMatlab__
#define __OutputToMatlab__

#include "DoubleArray1D.h"
#include "DoubleArray2D.h"

void outputToMatlab(DoubleArray2D& V, const char* filename);

#endif
 

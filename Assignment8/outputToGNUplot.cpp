//
//######################################################################
//
// outputToGNUplot.cpp : a routine that constructs of a gnuplot readable 
// data file that contains the the values of a scalar function of two 
// variables sampled on a cartesian product mesh, e.g. at the
// points (X(i),Y(j)) for i = 0 ... mPanels and j = 0 .. nPanels;
//                
// 
// Math 269B                                                 02/23/2010
//######################################################################
//
//

#include "DoubleArray1D.h"
#include "DoubleArray2D.h"
#include <math.h>

#include <iostream>
using namespace std;

void outputToGNUplot(DoubleArray1D& X, DoubleArray1D&Y, DoubleArray2D& F, const char* fileName)
{
	//
//  Open and then write to a file
//
    FILE* dataFile;

    if( (dataFile = fopen(fileName, "w+" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",fileName);
      return;
    }
 
	long mPanel = X.getSize()-1;
	long nPanel = Y.getSize()-1;

    long i; long j;
//
//  Output the data.
// 
    for(i = 0;  i <= mPanel;  i++)
    {
    for(j = 0;  j <= nPanel;  j++)
    {
    fprintf(dataFile,"%-10.5e %-10.5e %-10.5e \n",X(i),Y(j),F(i,j));
    }
    fprintf(dataFile,"\n");
    }
    
    fclose(dataFile);
} 

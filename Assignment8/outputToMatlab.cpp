#include "outputToMatlab.h"

//
//######################################################################
//
// outputToMatlab.cpp : a routine that constructs a Matlab readable 
// data file that contains the values of a scalar function of two 
// variables sampled on a Cartesian product mesh, e.g. at the
// grid points (i,j) for i = 0 ... mPanels and j = 0 .. nPanels;
//                
// 
// Math 269B                                                   02/23/10
//######################################################################
//
//

void outputToMatlab(DoubleArray2D& V, const char* filename)
{
//
//  Open and then write to a file
//
    FILE* dataFile;

    if( (dataFile = fopen(filename, "w+" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",filename);
      return;
    }

    long mPanel = V.getIndex1Size() - 1;
    long nPanel = V.getIndex2Size() - 1;

    long i; 
    long j;

    for(j = 0; j <= nPanel; j++)
    {
    for(i = 0; i <= mPanel; i++)
    {
    fprintf(dataFile,"%18.14e ",V(i,j));
    }
    fprintf(dataFile,"\n"); 
    }

    fclose(dataFile);
}
 

#ifndef __ProfileOutput__
#define __ProfileOutput__

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#include "DoubleArray1D.h"

class ProfileOutput
{
    public : 

    ProfileOutput(double xA, double xB)
    {
    this->xA = xA;
    this->xB = xB;
    };

 
    void output(DoubleArray1D& u, long outputIndex, char* fileNamePrefix, int outputFormat)
    {
    //
    // compose file name 
    //

    fileStringStream.str("");
    fileStringStream << fileNamePrefix;
    if     (outputIndex <= 9)  {fileStringStream << "00" << outputIndex << ".dat";}
    else if(outputIndex <= 99) {fileStringStream << "0" << outputIndex  << ".dat";}
    else                       {fileStringStream <<        outputIndex  << ".dat";}
//
//  Open and then write to a file
//
    FILE* dataFile;

    if( (dataFile = fopen((fileStringStream.str()).c_str(), "w+" )) == NULL )
    {
      printf( "The file %s could not be  opened\n",(fileStringStream.str()).c_str());
      return;
    }

    double mPanel = u.getSize() - 1;
    double dx     = (xB-xA)/double(mPanel);
  
    long i;
    double xPos;

    if(outputFormat != EXCEL) // MATLAB or GNUPLOT format space separated columns
    {
        for(i = 0; i < mPanel+1; i++)
        {
        xPos = xA + double(i)*dx;
        fprintf(dataFile,"%18.14e %18.14e",xPos,u(i));
        fprintf(dataFile," \n");
        }
    }
    else                     // EXCEL format (tab separated columns)
    {
        for(i = 0; i < mPanel+1; i++)
        {
        xPos = xA + double(i)*dx;
        fprintf(dataFile,"%18.14e \t %18.14e",xPos,u(i));
        fprintf(dataFile," \n");
        }
    }

    }

    double xA;                         // interval bounds 
    double xB; 

    enum {MATLAB,GNUPLOT,EXCEL};       // data file output types 
    ostringstream  fileStringStream;   // String stream creating output file titles

};
#endif
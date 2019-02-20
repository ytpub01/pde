#ifndef  __LaplaceOp2D__
#define  __LaplaceOp2D__
//
//######################################################
//                  Class LaplaceOp2D
//######################################################
//
// The apply(...) member function of this class applies the
// discrete 2D Laplace operator to interior values stored in a
// DoubleArray2D. The required Dirichlet boundary values
// are stored passed as additional arguments.
//
// Domain geometry
//
//   yB +-----------------+
//      |                 |
//      |                 | index j ^
//      |                 |
//   yA +-----------------+
//      xA               xB
//          index i >
//
//
//######################################################
// Initial Construction: Feb. 23, 2010
// 2/23/2010: Corrected indexing in "corner" apply
//######################################################
//
class LaplaceOp2D
{

public : 

//######################################################
// Constructor : Specify domain [xA,xB] x [yA, yB]
//######################################################

LaplaceOp2D(double xA, double xB, double yA, double yB)
{
    xMin = xA; xMax = xB;
    yMin = yA; yMax = yB;
}

//######################################################
//    void apply(DoubleArray2D& Vin,
//               DoubleArray1D& bdXA, DoubleArray1D& bdXB,
//               DoubleArray1D& bdYA, DoubleArray1D& bdYB,
//               DoubleArray2D& Vout);
//######################################################
//
// Applies the discrete 2D Laplace operator with boundary
// values specified in bdXA, bdXB, bdYA, bdYB
//
// On input  :  Vin, an array of size (mPanel-1, nPanel-1)
//              contains the values of a discrete function
//              associated with interior nodes of an equispaced
//              2D grid. The boundary conditions used for the 
//              the operator are the values specified in the arrays
//
//              bdXA = data at xA = array of size nPanel-1
//              bdXB = data at xB = array of size nPanel-1
//
//              bdYA = data at yA = array of size mPanel-1
//              bdYB = data at yB = array of size mPanel-1
//
//              Vout, and array of size (mPanel-1, nPanel-1)
//
// On return  : Vout, an array of size (mPanel-1, nPanel-1),
//              with  values set to the discrete Laplace
//              operator applied to the input values.
//
//              The contents of Vin are not changed by this routine.
// 
void apply(const DoubleArray2D& Vin,
		    const DoubleArray1D& bdXA, const DoubleArray1D& bdXB,
            const DoubleArray1D& bdYA, const DoubleArray1D& bdYB,
            DoubleArray2D& Vout)
{
    long mPanel = Vin.getIndex1Size() + 1;
    long nPanel = Vin.getIndex2Size() + 1;

    double hx = (xMax-xMin)/double(mPanel);
    double hy = (yMax-yMin)/double(nPanel);
    
    long i; long j;
    //
    // Array indices associated with interior data -
    //
    // i from 0 to mPanel-2
    // j from 0 to nPanel-2
    // 
    // Compute the discrete Laplacian at points not adjacent to a boundary
    // point.
    //

    for(j = 1; j < nPanel-2; j++)
    {
    for(i = 1; i < mPanel-2; i++)
    {
    Vout(i,j)  = (Vin(i+1,j) - 2.0*Vin(i,j) + Vin(i-1,j))/(hx*hx);
    Vout(i,j) += (Vin(i,j+1) - 2.0*Vin(i,j) + Vin(i,j-1))/(hy*hy);
    }}

    //
    // Evaluate discrete Laplace operator utilizing boundary data
    //

    //
    // Left and right edges at x =Corrected boundary indexing in Vplot packing 2/28/10 xA and x = xB
    //
    for(j = 1; j < nPanel-2; j++)
    {
    i = 0;
    Vout(i,j)  = (Vin(i+1,j) - 2.0*Vin(i,j) + bdXA(j))/(hx*hx);
    Vout(i,j) += (Vin(i,j+1) - 2.0*Vin(i,j) + Vin(i,j-1))/(hy*hy);

    i = mPanel-2;
    Vout(i,j)  = (bdXB(j) - 2.0*Vin(i,j) + Vin(i-1,j))/(hx*hx);
    Vout(i,j) += (Vin(i,j+1) - 2.0*Vin(i,j) + Vin(i,j-1))/(hy*hy);
    }


    //
    // Top and bottom edges at y = yA and y = yB
    //

    for(i = 1; i <  mPanel-2; i++)
    {
    j = 0;
    Vout(i,j)  = (Vin(i+1,j) - 2.0*Vin(i,j) + Vin(i-1,j))/(hx*hx);
    Vout(i,j) += (Vin(i,j+1) - 2.0*Vin(i,j) + bdYA(i))/(hy*hy);

    j = nPanel-2;
    Vout(i,j)  = (Vin(i+1,j) - 2.0*Vin(i,j) + Vin(i-1,j))/(hx*hx);
    Vout(i,j) += (bdYB(i) - 2.0*Vin(i,j) + Vin(i,j-1))/(hy*hy);
    }
    //
    // Corner Points
    //

	i = 0;
	j = 0;
	Vout(i,j)  = (Vin(i+1,j) - 2.0*Vin(i,j) + bdXA(j))/(hx*hx);
	Vout(i,j) += (Vin(i,j+1) - 2.0*Vin(i,j) + bdYA(i))/(hy*hy);

    i = mPanel-2;
	j = 0;
	Vout(i,j)  = (bdXB(j) - 2.0*Vin(i,j) + Vin(i-1,j))/(hx*hx);
	Vout(i,j) += (Vin(i,j+1) - 2.0*Vin(i,j) + bdYA(i))/(hy*hy);

	i = 0;
	j = nPanel-2;
	Vout(i,j)  = (Vin(i+1,j) - 2.0*Vin(i,j) + bdXA(j))/(hx*hx);
	Vout(i,j) += (bdYB(i) - 2.0*Vin(i,j) + Vin(i,j-1))/(hy*hy);

    i = mPanel-2;
    j = nPanel-2;
    Vout(i,j)  = (bdXB(j) - 2.0*Vin(i,j) + Vin(i-1,j))/(hx*hx);
    Vout(i,j) += (bdYB(i) - 2.0*Vin(i,j) + Vin(i,j-1))/(hy*hy);
}




    double xMax;  // Limits of the domain 
    double xMin;
    double yMax;
    double yMin; 
};

#endif
 

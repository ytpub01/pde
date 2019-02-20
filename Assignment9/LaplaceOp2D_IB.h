#ifndef  __LaplaceOp2D__
#define  __LaplaceOp2D__
//
//######################################################
//                 Class LaplaceOp2D_IB
//######################################################
//
// The apply(...) member function of this class applies the
// discrete 2D Laplace operator to values stored in a 
// DoubleArray2D that includes the Dirichlet boundary
// conditions.
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
//######################################################
// Feb. 23, 2010
//######################################################
//
class LaplaceOp2D_IB
{

public : 

//######################################################
// Constructor : Specify domain [xA,xB] x [yA, yB]
//######################################################

LaplaceOp2D_IB(double xA, double xB, double yA, double yB)
{
    xMin = xA; xMax = xB;
    yMin = yA; yMax = yB;
}

//######################################################
//    void apply(DoubleArray2D& Vin, DoubleArray2D& Vout)
//######################################################
//
// Applies the discrete 2D Laplace operator with Dirichlet 
// boundary conditions.
//
// On input  :  Vin, an array of size (mPanel+1, nPanel+1)
//              contains the values of a discrete function
//              associated with all nodes of an equispaced
//              2D grid. The boundary conditions used for the 
//              the operator are the values specified in the edge 
//              values of Vin. 
//
//              Vout, and array of size (mPanel+1, nPanel+1)
//
// On return  : Vout, an array of size (mPanel+1, nPanel+1), 
//              with interior values set to the discrete Laplace
//              operator applied to the input values, and 0 for the
//              boundary values. 
//// Domain geometry
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
//              The contents of Vin are not changed by this routine.
// 
void apply(const DoubleArray2D& Vin, DoubleArray2D& Vout)
{
    long mPanel = Vin.getIndex1Size() - 1;
    long nPanel = Vin.getIndex2Size() - 1;

    double hx = (xMax-xMin)/double(mPanel);
    double hy = (yMax-yMin)/double(nPanel);
    
    long i; long j;

    //
    // Array indices associated with data
    //
    // i from 0 to mPanel
    // j from 0 to nPanel
    //

    // 
    // compute the discrete Laplacian at interior points
    //

    for(j = 1; j <= nPanel-1; j++)
    {
    for(i = 1; i <= mPanel-1; i++)
    {
    Vout(i,j)  = (Vin(i+1,j) - 2.0*Vin(i,j) + Vin(i-1,j))/(hx*hx);
    Vout(i,j) += (Vin(i,j+1) - 2.0*Vin(i,j) + Vin(i,j-1))/(hy*hy);
    }}
    
    // Set edge values to zero

    for(j = 0; j <= nPanel; j++)   // left and right and bottom
    {
    Vout(0,j)       = 0.0;
    Vout(mPanel,j)  = 0.0;
    }

    for(i = 0; i <= mPanel; i++)   // top and bottom
    {
    Vout(i,0)       = 0.0;
    Vout(i,nPanel)  = 0.0;
    }
}

    double xMax;  // Limits of the domain 
    double xMin;
    double yMax;
    double yMin; 
};

#endif
 

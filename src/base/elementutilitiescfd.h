#ifndef incl_utilities_cfd_h
#define incl_utilities_cfd_h


#include "headersBasic.h"
#include "headersEigen.h"


inline void printMatrix(MatrixXd& AA)
{
    int ii, jj;
    printf("\n\n");
    for(ii=0;ii<AA.rows();ii++)
    {
        for(jj=0;jj<AA.cols();jj++)
           printf("\t%14.8f", AA(ii,jj));
        printf("\n");
    }
    printf("\n\n");

    return;
}


inline void printVector(VectorXd& AA)
{
    printf("\n\n");
    for(int ii=0;ii<AA.rows();ii++)
        printf("\t%6d\t%12.8f\n", ii, AA(ii));
    printf("\n\n");

   return;
}



inline void printVector(vector<int>&  vec)
{
    printf("\n\n");
    for(int ii=0;ii<vec.size();ii++)
        printf("\t%6d ", vec[ii]);
    printf("\n\n");

   return;
}


inline void printVector(double* data, int nn)
{
    printf("\n\n");
    for(int ii=0;ii<nn;ii++)
      printf("\t%6d\t%12.8f\n", ii, data[ii]);
    printf("\n\n");

   return;
}


void SetTimeParametersFluid(int tis, double rho, double dt, VectorXd& td);


int computeBasisFunctions2D(bool flag, int ETYPE, int degree, double* param, double* xNode, double* yNode, double*  N, double*  dN_dx, double* dN_dy, double&  Jac);


int computeBasisFunctions3D(bool flag, int ETYPE, int degree, double* param, double* xNode, double* yNode, double* zNode, double*  N, double*  dN_dx, double* dN_dy, double* dN_dz, double&  Jac);



#endif

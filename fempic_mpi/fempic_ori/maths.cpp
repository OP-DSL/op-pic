/*==============================================================================*
 * MATHS
 *------------------------------------------------------------------------------*
 * Maintainer: Ed Higgins <ed.higgins@york.ac.uk>
 * Based on `fem-pic.cpp` by Lubos Brieda 
 * See https://www.particleincell.com/2015/fem-pic/ for more information
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/


#include "maths.h"
#include <iostream>

/*computes determinant of a 4x4 matrix*/
double det4(double (*M)[4]) { 
    double M0[3][3];
    double M1[3][3];
    double M2[3][3];
    double M3[3][3];

    for (int i=0;i<3;i++) {
        M0[i][0]=M[i+1][1];
        M0[i][1]=M[i+1][2];
        M0[i][2]=M[i+1][3];

        M1[i][0]=M[i+1][0];
        M1[i][1]=M[i+1][2];
        M1[i][2]=M[i+1][3];

        M2[i][0]=M[i+1][0];
        M2[i][1]=M[i+1][1];
        M2[i][2]=M[i+1][3];

        M3[i][0]=M[i+1][0];
        M3[i][1]=M[i+1][1];
        M3[i][2]=M[i+1][2];
    }

    return M[0][0]*det3(M0) -
           M[0][1]*det3(M1) +
           M[0][2]*det3(M2) -
           M[0][3]*det3(M3);
}

/*computes determinant of a 3x3 matrix*/
double det3(double (*M)[3]) { 
    return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-
           M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+
           M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
}

/*helper functions for matrix math, y=A*x */
void matVecMultiply(double *y, double**A, double *x, int nu) { 
    // #pragma omp parallel for
    for (int i=0;i<nu;i++) {
        y[i] = 0;
        for (int j=0;j<nu;j++)
            y[i] += A[i][j]*x[j];
    }
}

/*computes y=v1-v2*/
void vecVecSubtract(double *y, double *v1, double *v2, int nu) {
    for (int i=0;i<nu;i++)
            y[i] = v1[i]-v2[i];
}

/*compute inverse of a 3x3 matrix using the adjugate method*/
void inverse(double M[3][3], double V[3][3]) 
{

    double a=M[0][0];
    double b=M[0][1];
    double c=M[0][2];
    double d=M[1][0];
    double e=M[1][1];
    double f=M[1][2];
    double g=M[2][0];
    double h=M[2][1];
    double i=M[2][2];

    V[0][0]=(e*i-f*h);
    V[1][0]=-(d*i-f*g);
    V[2][0]=(d*h-e*g);
    V[0][1]=-(b*i-c*h);
    V[1][1]=(a*i-c*g);
    V[2][1]=-(a*h-b*g);
    V[0][2]=(b*f-c*e);
    V[1][2]=-(a*f-c*d);
    V[2][2]=(a*e-b*d);
    double det = a*V[0][0]+b*V[1][0]+c*V[2][0];

    double Vmax = 0;
    for (int m=0;  m<3; m++) 
    {
        for (int n=0;  n<3; n++) 
        {
            Vmax = fabs(V[m][n]) > Vmax ? fabs(V[m][n]) : Vmax;
        }
    }

    double idet=0;
    if (fabs(Vmax) / fabs(det) > 1e12) 
    {
        std::cerr<<"Matrix is not invertible, |det M| = " << fabs(det) << 
            "! setting to [0]."<<std::endl;
    }
    else 
        idet=1/det;

    /*1/det*/
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            V[i][j]*=idet;
}

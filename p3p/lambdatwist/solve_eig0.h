#pragma once
#include "solve_cubic.h"
#include <Eigen/Dense>
#include <Eigen/Core>

/**
 * @brief eigwithknown0 eigen decomp of a matrix which has a 0 eigen value
 * @param x
 * @param E eigenvectors
 * @param L eigenvalues
 */
template<class T>
void eigwithknown0(Eigen::Matrix<T,3,3> x,
                   Eigen::Matrix<T,3,3>& E,
                   Eigen::Matrix<T,3,1>& L){
    // one eigenvalue is known to be 0.

    //the known one...
    L(2)=0;

    Eigen::Matrix<T,3,1>  v3(x(3)*x(7)- x(6)*x(4),
                   x(6)*x(1)- x(7)*x(0),
                   x(4)*x(0)- x(3)*x(1));
    v3.normalize();


    T x01_squared=x(0,1)*x(0,1);
    // get the two other...
    T b=- x(0,0) - x(1,1) - x(2,2);
    T c=- x01_squared - x(0,2)*x(0,2) - x(1,2)*x(1,2) +
            x(0,0)*(x(1,1) + x(2,2)) + x(1,1)*x(2,2);
    T e1,e2;
    root2real(b,c,e1,e2);

    if(std::abs(e1)<std::abs(e2))
        std::swap(e1,e2);
    L(0)=e1;
    L(1)=e2;


    T mx0011=-x(0,0)*x(1,1);
    T prec_0 = x(0,1)*x(1,2) - x(0,2)*x(1,1);
    T prec_1 = x(0,1)*x(0,2) - x(0,0)*x(1,2);


    T e=e1;
    T tmp=1.0/(e*(x(0,0) + x(1,1)) + mx0011 - e*e + x01_squared);
    T a1= -(e*x(0,2) + prec_0)*tmp;
    T a2= -(e*x(1,2) + prec_1)*tmp;
    T rnorm=((T)1.0)/std::sqrt(a1*a1 +a2*a2 + 1.0);
    a1*=rnorm;
    a2*=rnorm;
    Eigen::Matrix<T,3,1>  v1(a1,a2,rnorm);


    T tmp2=1.0/(e2*(x(0,0) + x(1,1)) + mx0011 - e2*e2 + x01_squared);
    T a21= -(e2*x(0,2) + prec_0)*tmp2;
    T a22= -(e2*x(1,2) + prec_1)*tmp2;
    T rnorm2=1.0/std::sqrt(a21*a21 +a22*a22 +1.0);
    a21*=rnorm2;
    a22*=rnorm2;
    Eigen::Matrix<T,3,1>  v2(a21,a22,rnorm2);


E<<v1[0],v2[0],v3[0],  v1[1],v2[1],v3[1],  v1[2],v2[2],v3[2];

}
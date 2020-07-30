#pragma once
/* 
This code is directly taken from the this paper: 
https://openaccess.thecvf.com/content_ECCV_2018/papers/Mikael_Persson_Lambda_Twist_An_ECCV_2018_paper.pdf
The link to the paper also contains a link to the git repo containing original source code
*/
#include "solve_cubic.h"
#include "solve_eig0.h"
#include "refine_lambda.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

using std::endl;using std::cout;

template<class T, int refinement_iterations=5>
int p3p_lambdatwist(const Eigen::Matrix<T, 3, 1>& world_pt1,
                    const Eigen::Matrix<T, 3, 1>& world_pt2,
                    const Eigen::Matrix<T, 3, 1>& world_pt3,
                    const Eigen::Matrix<T, 3, 1>& cam_pt1,
                    const Eigen::Matrix<T, 3, 1>& cam_pt2,
                    const Eigen::Matrix<T, 3, 1>& cam_pt3,
                    Eigen::Matrix<Eigen::Matrix<T,3,3>,4,1>& Rs,
                    Eigen::Matrix<Eigen::Matrix<T, 3, 1>,4,1>& Ts){

    
    Eigen::Matrix<T, 3, 1> y1 = world_pt1;
    Eigen::Matrix<T, 3, 1> y2 = world_pt2;
    Eigen::Matrix<T, 3, 1> y3 = world_pt3;
    
    // normalize the length of ys, we could expect it, but lets not...
    y1.normalize();
    y2.normalize();
    y3.normalize();


    T b12=-2.0*(y1.dot(y2));
    T b13=-2.0*(y1.dot(y3));
    T b23=-2.0*(y2.dot(y3));


    // implicit creation of Eigen::Matrix<T, 3, 1> can be removed
    Eigen::Matrix<T, 3, 1> d12=cam_pt1-cam_pt2;
    Eigen::Matrix<T, 3, 1> d13=cam_pt1-cam_pt3;
    Eigen::Matrix<T, 3, 1> d23=cam_pt2-cam_pt3;
    Eigen::Matrix<T, 3, 1> d12xd13(d12.cross(d13));


    T a12=d12.norm();
    T a13=d13.norm();
    T a23=d23.norm();


    T c31=-0.5*b13;
    T c23=-0.5*b23;
    T c12=-0.5*b12;
    T blob=(c12*c23*c31 - 1.0);

    T s31_squared=1.0-c31*c31;
    T s23_squared=1.0-c23*c23;
    T s12_squared=1.0-c12*c12;


    T p3 = (a13*(a23*s31_squared - a13*s23_squared));

    T p2 = 2.0*blob*a23*a13 + a13*(2.0*a12 + a13)*s23_squared + a23*(a23 -     a12)*s31_squared;

    T p1 = a23*(a13 - a23)*s12_squared - a12*a12*s23_squared - 2.0*a12*(blob*a23 +    a13*s23_squared);

    T p0 = a12*(a12*s23_squared - a23*s12_squared);

    T g=0;

    {
        p3=1.0/p3;
        p2*=p3;
        p1*=p3;
        p0*=p3;

        g=cubick<T>(p2,p1,p0);
    }

    T A00=a23*(1.0- g);
    T A01=(a23*b12)*0.5;
    T A02=(a23*b13*g)*(-0.5);
    T A11=a23 - a12 + a13*g;
    T A12=b23*(a13*g - a12)*0.5;
    T A22=g*(a13 - a23) - a12;

    Eigen::Matrix<T,3,3> A;
    A << A00,A01,A02,
          A01,A11,A12,
          A02,A12,A22;

    // get sorted eigenvalues and eigenvectors given that one should be zero...
    Eigen::Matrix<T,3,3> V;
    Eigen::Matrix<T,3,1> L;
    eigwithknown0<T>(A,V,L);

    T v=std::sqrt(std::max(T(0),-L(1)/L(0)));


    int valid=0;
    Eigen::Matrix<Eigen::Matrix<T,3,1>,4,1> Ls;
    

    {
        T s=v;

        T w2=T(1.0)/( s*V(1) - V(0));
        T w0=(V(3) - s*V(4))*w2;
        T w1=(V(6) - s*V(7))*w2;



        T a=T(1.0)/((a13 - a12)*w1*w1 - a12*b13*w1 - a12);
        T b=(a13*b12*w1 - a12*b13*w0 - T(2.0)*w0*w1*(a12 - a13))*a;
        T c=((a13 - a12)*w0*w0 + a13*b12*w0 + a13)*a;



        if(b*b -4.0*c>=0 ){
            T tau1,tau2;
            root2real<T>(b,c,tau1,tau2);
            if(tau1>0){
                T tau=tau1;
                T d=a23/(tau*(b23 + tau) + T(1.0));

                T l2=std::sqrt(d);
                T l3=tau*l2;

                T l1=w0*l2 +w1*l3;
                if(l1>=0){

                    Ls[valid]={l1,l2,l3};

                    ++valid;
                }

            }
            if(tau2>0){
                T tau=tau2;
                T d=a23/(tau*(b23 + tau) + T(1.0));

                T l2=std::sqrt(d);
                T l3=tau*l2;
                T l1=w0*l2 +w1*l3;
                if(l1>=0){
                    Ls[valid]={l1,l2,l3};
                    ++valid;
                }

            }
        }
    }

    {
        T s=-v;
        T w2=T(1.0)/( s*V(0,1) - V(0,0));
        T w0=(V(1,0) - s*V(1,1))*w2;
        T w1=(V(2,0) - s*V(2,1))*w2;

        T a=T(1.0)/((a13 - a12)*w1*w1 - a12*b13*w1 - a12);
        T b=(a13*b12*w1 - a12*b13*w0 - T(2.0)*w0*w1*(a12 - a13))*a;
        T c=((a13 - a12)*w0*w0 + a13*b12*w0 + a13)*a;


        if(b*b -4.0*c>=0){
            T tau1,tau2;

            root2real<T>(b,c,tau1,tau2);
            if(tau1>0) {
                T tau=tau1;
                T d=a23/(tau*(b23 + tau) + T(1.0));
                if(d>0){
                  T l2=std::sqrt(d);

                  T l3=tau*l2;

                  T l1=w0*l2 +w1*l3;
                  if(l1>=0){
                      Ls[valid]={l1,l2,l3};
                      ++valid;
                  }
                }
            }
            if(tau2>0){
                T tau=tau2;
                T d=a23/(tau*(b23 + tau) + T(1.0));
                if(d>0){
                  T l2=std::sqrt(d);

                  T l3=tau*l2;

                  T l1=w0*l2 +w1*l3;
                  if(l1>=0){
                      Ls[valid]={l1,l2,l3};
                      ++valid;
                  }
                }
            }
        }
    }

    for(int i=0;i<valid;++i){              
        gauss_newton_refineL<T,refinement_iterations>(Ls[i],a12,a13,a23,b12,b13,b23);        
    }

    Eigen::Matrix<T, 3, 1> ry1,ry2,ry3;
    Eigen::Matrix<T, 3, 1> yd1;
    Eigen::Matrix<T, 3, 1> yd2;
    Eigen::Matrix<T, 3, 1> yd1xd2;
    Eigen::Matrix<T,3,3> X;
    X << d12(0),d13(0),d12xd13(0),
         d12(1),d13(1),d12xd13(1),
         d12(2),d13(2),d12xd13(2);
    X=X.inverse().eval();

    for(int i=0;i<valid;++i){

        // compute the rotation:
        ry1=y1*Ls(i)(0);
        ry2=y2*Ls(i)(1);
        ry3=y3*Ls(i)(2);

        yd1=ry1-ry2;
        yd2=ry1-ry3;
        yd1xd2=yd1.cross(yd2);

        Eigen::Matrix<T,3,3> Y;
        Y << yd1(0),yd2(0),yd1xd2(0),
              yd1(1),yd2(1),yd1xd2(1),
              yd1(2),yd2(2),yd1xd2(2);


        Rs[i]=Y*X; // probably not needed cause data can be retrieved other ways

        Ts[i]=(ry1 - Rs[i]*cam_pt1 ); // probably not needed cause data can be retrieved other ways
    
    }

    return valid;


}


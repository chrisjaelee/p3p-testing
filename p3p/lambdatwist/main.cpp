#include "test_pnpransac.h"
// #include "lambdatwist_pnpransac.h"
// #include "utils/cvl/matrix.h"
// #include <iostream>
// #include <Eigen/Dense>
// #include <opencv2/opencv.hpp>

int main(){

     test_pnpransac();
    // Eigen::Matrix3d image_pts;
    // image_pts << 0.82561,  0.258205,  1,
    //              0.811942, -0.455175, 1,
    //              0.922062, 0.421522,  1;
    // image_pts.transposeInPlace();
    
    // Eigen::Matrix<double,3,3>  world_pts;
    // world_pts << 12.8451, 4.34078, 24.2388,
    //              3.6176,  32.5014, 1.97276,
    //              34.6513, 36.495,  39.9275;

    // // Eigen::Matrix<double,4,3>  worldPtsForReprojection;
    // // worldPtsForReprojection << 12.8451, 4.34078, 24.2388,
    // //                            3.6176,  32.5014, 1.97276,
    // //                            34.6513, 36.495,  39.9275,
    // //                            1,       1,       1;
    // // worldPtsForReprojection.block(0,0,3,3) = (worldPtsForReprojection.block(0,0,3,3).transpose()).eval();

    // // image_pts = (worldPtsForReprojection.block(0,0,3,3).transpose()).eval();
                 
    // Eigen::Matrix<double,3,3> camFromWorldRotation;
    // camFromWorldRotation << -0.867744, -0.472055,  0.155513,
    //                         -0.494428,  0.851753, -0.173375,
    //                         -0.050616, -0.227335,   -0.9725;

    // Eigen::Matrix<double,3,1> camFromWorldTranslation;
    // camFromWorldTranslation << -21.3241,
    //                            0.275653,
    //                            0.291536;


    // Eigen::Matrix<double,4,4> camFromWorldTransform = Eigen::MatrixXd::Identity(4,4); 
    // camFromWorldTransform.block(0,0,3,3) = camFromWorldRotation;
    // camFromWorldTransform.block(0,3,3,1) = camFromWorldTranslation;

    // Eigen::Matrix<double,4,3> reprojectedPts = camFromWorldTransform*worldPtsForReprojection;
    // for(int i=0;i<3;i++){
    //     reprojectedPts.col(i) = (reprojectedPts.col(i)/reprojectedPts(2,i));
    // }
    // std::cout << reprojectedPts;

    // cvl::Vector<cvl::Matrix<double,3,3>,4> Rs;
    // cvl::Vector<cvl::Vector<double,3>,4> Ts;

    // Eigen::Matrix<double,3,1> temp;

    // cvl::Matrix<double, 3, 1> x0;
    // cvl::Matrix<double, 3, 1> x1;
    // cvl::Matrix<double, 3, 1> x2;
    // temp = image_pts.col(0);
    // x0 = eigenToCVL<double,3,1>(temp);
    // std::cout << "\nImage Points:  " << x0.transpose(); 
    // temp = image_pts.col(1);
    // x1 = eigenToCVL<double,3,1>(temp);
    // std::cout << "Image Points:  " << x1.transpose(); 
    // temp = image_pts.col(2);
    // x2 = eigenToCVL<double,3,1>(temp);
    // std::cout << "Image Points:  " << x2.transpose(); 

    // cvl::Matrix<double, 3, 1> y0;
    // cvl::Matrix<double, 3, 1> y1;
    // cvl::Matrix<double, 3, 1> y2;
    // temp = world_pts.col(0);
    // y0 = eigenToCVL<double,3,1>(temp);
    // std::cout << "World Points:  " << y0.transpose(); 
    // temp = world_pts.col(1);
    // y1 = eigenToCVL<double,3,1>(temp);
    // std::cout << "World Points:  " << y1.transpose(); 
    // temp = world_pts.col(2);
    // y2 = eigenToCVL<double,3,1>(temp);
    // std::cout << "World Points:  " << y2.transpose(); 

    // int valid = cvl::p3p_lambdatwist(x0,x1,x2,y0,y1,y2,
    //                                  Rs, Ts);

    // std::cout << "LAMBDATWIST\n";
    // for(int i=0;i<valid;i++){
    //     std::cout << "ROTATION: \n" << to_rodrigues(cvlToEigen(Rs[i]));
    //     std::cout << "\nTRANSLATION: \n" << Ts[i] << std::endl;
    // }
    
    // //CV STUFF
    // std::vector<cv::Point2d> cvImagePoints;
    // cvImagePoints.push_back(cv::Point2d(x0[0],x0[1]));
    // cvImagePoints.push_back(cv::Point2d(x1[0],x1[1]));
    // cvImagePoints.push_back(cv::Point2d(x2[0],x2[1]));

    // std::vector<cv::Point3d> cvWorldPoints;
    // cvWorldPoints.push_back(cv::Point3d(y0[0],y0[1],y0[2]));
    // cvWorldPoints.push_back(cv::Point3d(y1[0],y1[1],y1[2]));
    // cvWorldPoints.push_back(cv::Point3d(y2[0],y2[1],y2[2]));
    
    // std::vector<cv::Mat> Rs_cv; 
    // std::vector<cv::Mat> Ts_cv;

    // cv::Mat camera_matrix = (cv::Mat_<double>(3,3) << 1, 0, 0, 
    //                                                   0, 1, 0,
    //                                                   0, 0, 1);

    // cv::solveP3P(cvWorldPoints, cvImagePoints, camera_matrix, cv::noArray(), 
    //              Rs_cv, Ts_cv, cv::SOLVEPNP_P3P);
    
    // std::cout << "LAMBDATWIST\n";
    // for(int i=0;i<Rs_cv.size();i++){
    //     std::cout << "ROTATION: \n" << Rs_cv[i];
    //     std::cout << "\nTRANSLATION: \n" << Ts_cv[i] << std::endl;
    // }
    
    

// gt rotation:
// -0.76421
// 2.91931
// -0.316846
// gt translation:
// -21.3241
// 0.275653
// 0.291536


    return 0; 
} 
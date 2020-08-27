#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>
#include "generator.h"
#include "utils/cvl/matrix.h"
#include "lambdatwist_pnpransac.h"

#define pi 3.14159265
#define big_number 99999999


Eigen::Matrix<double,3,1> cvVecToEigen(cv::Vec3d values){
    Eigen::Matrix<double,3,1> output;
    for(int i=0;i<3;i++){
        output[i] = values[i];
    }
    return output;
}

double findError(Eigen::Matrix<double,3,1> estimate, Eigen::Matrix<double,3,1> actual){
    double error = 0;
    Eigen::Matrix<double,3,1> difs = estimate-actual;
    for(int i=0;i<3;i++){
        error+=abs(difs[i]);
    }

    return error/2;
}

void test_pnpransac(const int& tests=1000, const bool& verbose=false){
    srand(1234);

    Eigen::Matrix<double,3,3> k;
    k << 1016, 0,   933,
         0,    1031, 521,
         0,    0,   1;

    Eigen::MatrixXd pts_in_world;

    Eigen::MatrixXd pts_in_cam;

    Eigen::Vector4d cam_in_world;
    
    Eigen::Matrix<double,4,4> cam_from_world_transform;

    Eigen::Matrix3d temp;

    Eigen::Matrix<double,3,1> gt_rotation;
    Eigen::MatrixXd gt_translation;

    std::vector<cv::Point2d> image_pts;
    std::vector<cv::Point3d> world_pts;

    cv::Vec3d Rs; 
    cv::Vec3d Ts;
    cv::Vec3d r;
    cv::Vec3d t;

    cv::Mat camera_matrix = (cv::Mat_<double>(3,3) << k(0,0), 0, k(0,2), 
                                                      0, k(1,1), k(1,2),
                                                      0, 0, 1);

    double twist_time = 0;
    double twist_rError = 0;
    double twist_tError = 0;

    double cv_time = 0;
    double cv_rError = 0;
    double cv_tError = 0;
    
    for(int i=0;i<tests;i++){
        while(true){
            pts_in_world = gen_pt_cloud<3000>();

            cam_from_world_transform = world_pts_in_cam(pts_in_world, pts_in_cam, k);
        
            for(int j=0;j<pts_in_cam.cols();j++){
                image_pts.push_back(cv::Point2d(pts_in_cam(0,j), pts_in_cam(1,j)));
                world_pts.push_back(cv::Point3d(pts_in_world(0,j), pts_in_world(1,j), pts_in_world(2,j)));
            }

            if(image_pts.size()>50){
                break;
            }
        }
        
        temp = cam_from_world_transform.block(0,0,3,3);
        Eigen::Matrix<double,3,1> gt_rotation = to_rodrigues(temp);
        Eigen::Matrix<double,3,1> gt_translation = cam_from_world_transform.block(0,3,3,1);

        auto t1 = std::chrono::high_resolution_clock::now();
        twistPnPRansac(world_pts,image_pts,camera_matrix,Rs,Ts);
        auto t2 = std::chrono::high_resolution_clock::now();
        twist_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
        twist_time = twist_time/2;

        twist_rError = (twist_rError + findError(cvVecToEigen(Rs),gt_rotation))/2;
        twist_tError = (twist_tError + findError(cvVecToEigen(Ts),gt_translation))/2;


        auto t3 = std::chrono::high_resolution_clock::now();
        cv::solvePnPRansac(world_pts,image_pts,camera_matrix,cv::noArray(),r,t,2);
        auto t4 = std::chrono::high_resolution_clock::now();
        cv_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t4 - t3 ).count();
        cv_time = cv_time/2;

        cv_rError = (cv_rError + findError(cvVecToEigen(r),gt_rotation))/2;
        cv_tError = (cv_tError + findError(cvVecToEigen(t),gt_translation))/2;

        if(verbose){
            std::cout << "\n\n------------\nTEST: " << i+1 << "/" << tests;
            std::cout << "\n\ntwist rError:\n" << twist_rError;
            std::cout << "\ntwist tError:\n" << twist_tError;
            std::cout << "\ntwist timer:\n" << twist_time;

            std::cout << "\n\ncv rError:\n" << cv_rError;
            std::cout << "\ncv tError:\n" << cv_tError;
            std::cout << "\ncv timer:\n" << cv_time;
        }
    }

    // std::cout << "\n------------\ntwist rError:\n" << twist_rError;
    // std::cout << "\ntwist tError:\n" << twist_tError;
    // std::cout << "\ntwist timer:\n" << twist_time;

    // std::cout << "\n\ncv rError:\n" << cv_rError;
    // std::cout << "\ncv tError:\n" << cv_tError;
    // std::cout << "\ncv timer:\n" << cv_time;
}
#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>
// #include "lambdatwist.p3p.h"
#include "generator.h"
#include "utils/cvl/matrix.h"
#include "lambdatwist_pnpransac.h"

#define pi 3.14159265
#define big_number 99999999




void test_pnpransac(const int& tests=1000, const bool& verbose=false){

    Eigen::Matrix<double,3,3> k;
    k << 1016, 0,   933,
         0,    1031, 521,
         0,    0,   1;

    Eigen::MatrixXd pts_in_world;

    Eigen::MatrixXd pts_in_cam;

    Eigen::Vector4d cam_in_world;
    
    Eigen::Matrix<double,4,4> transform;

    Eigen::Matrix3d temp;

    Eigen::Matrix<double,3,1> gt_rotation;
    Eigen::MatrixXd gt_translation;

    std::vector<cv::Point2d> image_pts;
    std::vector<cv::Point3d> world_pts;

    std::vector<cv::Mat> Rs; 
    std::vector<cv::Mat> Ts;
    cv::Vec3d r;
    cv::Vec3d t;

    cv::Mat camera_matrix = (cv::Mat_<double>(3,3) << k(0,0), 0, k(0,2), 
                                                      0, k(1,1), k(1,2),
                                                      0, 0, 1);
    // for(int i=0;i<tests;i++){
        while(true){
            pts_in_world = gen_pt_cloud<3000>();

            transform = world_pts_in_cam(pts_in_world, pts_in_cam, cam_in_world, k);
        
            for(int j=0;j<pts_in_cam.cols();j++){
                image_pts.push_back(cv::Point2d(pts_in_cam(0,j), pts_in_cam(1,j)));
                world_pts.push_back(cv::Point3d(pts_in_world(0,j), pts_in_world(1,j), pts_in_world(2,j)));
            }

            if(image_pts.size()>50){
                break;
            }
        }
        temp = transform.block(0,0,3,3);
        std::cout << to_rodrigues(temp) << "\n\n";
        

        // temp = transform.block(0,0,3,3);
        // gt_rotation = to_rodrigues(temp);
        // gt_translation = transform.block(0,3,3,1);

        twistPnPRansac(world_pts,image_pts,camera_matrix,Rs,Ts);
        cv::solvePnPRansac(world_pts,image_pts,camera_matrix,cv::noArray(),r,t);
        std::cout << "\n\ncv output" << r;
    // }
}
#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include "lambdatwist.p3p.h"
#include "generator.h"
#include "utils/cvl/matrix.h"

#define pi 3.14159265
 
using namespace Eigen;

template <class T>
Eigen::Matrix<T, 3, 1> to_rodrigues(const Eigen::Matrix<T, 3, 3>& data){
    Eigen::AngleAxis<double> aa (data);
    return aa.axis() * aa.angle();
}


template <class T, unsigned int rows, unsigned int cols>
T find_error(const Eigen::Matrix<T,rows,cols>& estimate, const Eigen::Matrix<T,rows,cols>& actual){
    T error{0};
    Eigen::Matrix<T,rows,cols> temp = estimate - actual;
    for(int i=0;i<rows;i++){
        for(int j=0;i<cols;i++){
            error += (std::abs(temp(i,j))/2);
        }
    }
    return error;
}


/*--------------------------
--------LAMBDA TWIST--------
--------------------------*/
//CVL::Matrix TO EIGEN
template <class T, unsigned int rows, unsigned int cols>
Eigen::Matrix<T, rows, cols> cvl_to_eigen(const cvl::Matrix<T, rows, cols>& data){
    Eigen::Matrix<T, rows, cols> eigen_matrix;
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            eigen_matrix(i,j) = data(i,j);
        }
    }
    return eigen_matrix;
}

Eigen::Matrix<double,3,1> run_lambdatwist(const Eigen::Matrix<double,4,3>& pts_3d,
                     const Eigen::MatrixXd& pts_2d,
                     const Eigen::MatrixXd& gt_rotation,
                     const Eigen::MatrixXd& gt_translation){

    Eigen::Matrix<double,3,1> output = Eigen::MatrixXd::Zero(3,1);

    cvl::Matrix<double, 3, 1> x0_cvl;
    cvl::Matrix<double, 3, 1> x1_cvl;
    cvl::Matrix<double, 3, 1> x2_cvl;
    x0_cvl = {pts_3d(0,0), pts_3d(1,0), pts_3d(2,0)};
    x1_cvl = {pts_3d(0,1), pts_3d(1,1), pts_3d(2,1)};
    x2_cvl = {pts_3d(0,2), pts_3d(1,2), pts_3d(2,2)};

    cvl::Matrix<double, 3, 1> y0_cvl;
    cvl::Matrix<double, 3, 1> y1_cvl;
    cvl::Matrix<double, 3, 1> y2_cvl;
    y0_cvl = {pts_2d(0,0), pts_2d(1,0), 1};
    y1_cvl = {pts_2d(0,1), pts_2d(1,1), 1};
    y2_cvl = {pts_2d(0,2), pts_2d(1,2), 1};

    cvl::Vector<cvl::Matrix<double,3,3>,4> Rs;
    cvl::Vector<cvl::Vector<double,3>,4> Ts;

    auto t1 = std::chrono::high_resolution_clock::now();
    int valid = cvl::p3p_lambdatwist(y0_cvl, y1_cvl, y2_cvl, x0_cvl, x1_cvl, x2_cvl, Rs, Ts);
    auto t2 = std::chrono::high_resolution_clock::now();

    output(2) = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

    // finding error!
    output(0) = find_error<double,3,1>(to_rodrigues(cvl_to_eigen(Rs(0))), gt_rotation);
    output(1) = find_error<double,3,1>(cvl_to_eigen(Ts(0)), gt_translation);

    double temp_r_error{0};
    double temp_t_error{0};

    for(int i = 1; i<valid; i++){
        temp_r_error = find_error<double,3,1>(to_rodrigues(cvl_to_eigen(Rs(i))), gt_rotation);
        temp_t_error = find_error<double,3,1>(cvl_to_eigen(Ts(i)), gt_translation);

        if(output(0)+output(1) > temp_r_error+temp_t_error){
            output(0) = temp_r_error;
            output(1) = temp_t_error;
        }
    }

    return output;
}


/*--------------------------
--------OPENCV P3P----------
--------------------------*/

template <class T, unsigned int rows, unsigned int cols>
Eigen::Matrix<T, rows, cols> cv_to_eigen(const cv::Mat& data){
    Eigen::Matrix<T, rows, cols> eigen_matrix;
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            eigen_matrix(i,j) = data.at<T>(i,j);
        }
    }
    return eigen_matrix;
}

Eigen::Matrix<double,3,1> run_cv_p3p(const Eigen::Matrix<double,4,3>& pts_3d,
                                     const Eigen::MatrixXd& pts_2d,
                                     Eigen::Matrix<double,3,1> gt_rotation,
                                     Eigen::Matrix<double,3,1> gt_translation,
                                     Eigen::Matrix<double,3,3> k){

    Eigen::Matrix<double,3,1> output = Eigen::MatrixXd::Zero(3,1);

    Eigen::Matrix<double, 3, 1> y0;
    Eigen::Matrix<double, 3, 1> y1;
    Eigen::Matrix<double, 3, 1> y2;
    y0 = {pts_3d(0,0), pts_3d(1,0), pts_3d(2,0)};
    y1 = {pts_3d(0,1), pts_3d(1,1), pts_3d(2,1)};
    y2 = {pts_3d(0,2), pts_3d(1,2), pts_3d(2,2)};

    Eigen::Matrix<double, 3, 1> x0;
    Eigen::Matrix<double, 3, 1> x1;
    Eigen::Matrix<double, 3, 1> x2;
    x0 = {pts_2d(0,0), pts_2d(1,0), 1};
    x1 = {pts_2d(0,1), pts_2d(1,1), 1};
    x2 = {pts_2d(0,2), pts_2d(1,2), 1};

    std::vector<cv::Point2d> image_points;

    Eigen::Matrix<double,3,1> x0_pixel = k*x0;
    Eigen::Matrix<double,3,1> x1_pixel = k*x1;
    Eigen::Matrix<double,3,1> x2_pixel = k*x2;
    image_points.push_back(cv::Point2d(x0_pixel[0], x0_pixel[1]));    
    image_points.push_back(cv::Point2d(x1_pixel[0], x1_pixel[1]));    
    image_points.push_back(cv::Point2d(x2_pixel[0], x2_pixel[1]));    
    
    // 3D model points.
    std::vector<cv::Point3d> model_points;
    model_points.push_back(cv::Point3d(y0[0], y0[1], y0[2]));  
    model_points.push_back(cv::Point3d(y1[0], y1[1], y1[2]));  
    model_points.push_back(cv::Point3d(y2[0], y2[1], y2[2]));
    
    // Camera internals
    cv::Mat camera_matrix = (cv::Mat_<double>(3,3) << k(0,0), 0, k(0,2), 
                                                      0, k(1,1), k(1,2),
                                                      0, 0, 1);
    
    std::vector<cv::Mat> Rs; 
    std::vector<cv::Mat> Ts;
    
    auto t1 = std::chrono::high_resolution_clock::now();
    cv::solveP3P(model_points, image_points, camera_matrix, cv::noArray(), 
                 Rs, Ts, cv::SOLVEPNP_P3P);
    auto t2 = std::chrono::high_resolution_clock::now();

    output(2) = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

    // finding error!
    cv_to_eigen<double,3,1>(Rs[0]);
    // output(0) = find_error<double,3,1>(cv_to_eigen<double,3,1>(Rs[0]), gt_rotation);
    // output(1) = find_error<double,3,1>(cv_to_eigen<double,3,1>(Ts[0]), gt_translation);

    // double temp_r_error{0};
    // double temp_t_error{0};
    
    // for(int i = 1; i<Rs.size(); i++){
    //     temp_r_error = find_error<double,3,1>(cv_to_eigen<double,3,1>(Rs[i]), gt_rotation);
    //     temp_t_error = find_error<double,3,1>(cv_to_eigen<double,3,1>(Ts[i]), gt_translation);

    //     if(output(0)+output(1) > temp_r_error+temp_t_error){
    //         output(0) = temp_r_error;
    //         output(1) = temp_t_error;
    //     }
    // }
    // return output;
}



void test_p3p(int tests=1000, const bool& verbose=false){
    Eigen::Matrix<double,3,3> k;
    k << 1016, 0,   933,
         0,    1031, 521,
         0,    0,   1;

    Eigen::Matrix<double,4,3> triangle_in_world;
    triangle_in_world << 2, 7, 4.5,
                         1, 1, 5.33,
                         20, 20, 20,
                         1, 1, 1;

    Eigen::MatrixXd triangle_in_cam;

    Eigen::Vector4d cam_in_world;
    
    Eigen::Matrix<double,4,4> transform;

    Eigen::Matrix<double,3,1> lambda_error = Eigen::MatrixXd::Zero(3,1);
    Eigen::Matrix<double,3,1> cv_error = Eigen::MatrixXd::Zero(3,1);
    Eigen::Matrix3d temp;

    Eigen::Matrix<double,3,1> gt_rotation;
    Eigen::MatrixXd gt_translation;

    for(int i=0;i<tests;i++){
        transform = world_pts_in_cam(triangle_in_world, triangle_in_cam, cam_in_world, k);

        temp = transform.block(0,0,3,3);
        gt_rotation = to_rodrigues(temp);
        gt_translation = transform.block(0,3,3,1);
        
        lambda_error = lambda_error + run_lambdatwist(triangle_in_world, triangle_in_cam, gt_rotation, gt_translation);
        lambda_error(2) = lambda_error(2)/2;
        cv_error = cv_error + run_cv_p3p(triangle_in_world, triangle_in_cam, gt_rotation, gt_translation, k);
        cv_error(2) = cv_error(2)/2;

    }
    if(verbose){
        std::cout << std::endl << "---LAMBDATWIST---" << std::endl;
        std::cout << lambda_error;
        std::cout << std::endl << "-----OPEN CV-----" << std::endl;
        std::cout << cv_error << std::endl;
    }
}
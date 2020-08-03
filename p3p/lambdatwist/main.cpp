#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>
#include "lambdatwist.p3p.h"
#include "generator.h"
#include "utils/cvl/matrix.h"
// #include "p3p_generator.h"

#define pi 3.14159265
 
using namespace Eigen;

template <class T, unsigned int rows, unsigned int cols>
Eigen::Matrix<T, rows, cols> to_eigen_matrix(const cvl::Matrix<T, rows, cols>& data){
    Eigen::Matrix<T, rows, cols> eigen_matrix;
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            eigen_matrix(i,j) = data(i,j);
        }
    }
    return eigen_matrix;
}

template <class T>
Eigen::Matrix<T, 3, 1> to_rodrigues(const Eigen::Matrix<T, 3, 3>& data){
    Eigen::AngleAxis<double> aa (data);
    return aa.axis() * aa.angle();
}

void test_p3p(){
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
    // cam_in_world << .84, 1.2, 1, 1;
    cam_in_world << 0,0,0,1;
    
    world_pts_in_cam(triangle_in_world, triangle_in_cam, cam_in_world, k, true);
    std::cout << std::endl;
    std::cout << "Triangle in cam: " << std::endl;
    std::cout << triangle_in_cam << std::endl;
    std::cout << std::endl;

     /*-------------------------
    --------OPENCV PNP----------
    --------------------------*/

    Eigen::Matrix<double, 3, 1> y0;
    Eigen::Matrix<double, 3, 1> y1;
    Eigen::Matrix<double, 3, 1> y2;
    y0 = {triangle_in_world(0,0), triangle_in_world(1,0), triangle_in_world(2,0)};
    y1 = {triangle_in_world(0,1), triangle_in_world(1,1), triangle_in_world(2,1)};
    y2 = {triangle_in_world(0,2), triangle_in_world(1,2), triangle_in_world(2,2)};

    Eigen::Matrix<double, 3, 1> x0;
    Eigen::Matrix<double, 3, 1> x1;
    Eigen::Matrix<double, 3, 1> x2;
    x0 = {triangle_in_cam(0,0), triangle_in_cam(1,0), 1};
    x1 = {triangle_in_cam(0,1), triangle_in_cam(1,1), 1};
    x2 = {triangle_in_cam(0,2), triangle_in_cam(1,2), 1};

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

    std::cout << "Camera Matrix " << std::endl << camera_matrix << std::endl ;
    
    std::vector<cv::Mat> rotation_vector; 
    std::vector<cv::Mat> translation_vector;

    cv::solveP3P(model_points, image_points, camera_matrix, cv::noArray(), 
                 rotation_vector, translation_vector, cv::SOLVEPNP_P3P);
    // double duration_cv = funcTime(cv::solvePnPRansac, model_points, image_points, 
    //                               camera_matrix, dist_coeffs, rotation_vector, translation_vector,
    //                               false, CV_P3P);
    std::cout << std::endl << "---OPENCV---" << std::endl;
    // std::cout << std::endl << duration_cv << " nanoseconds" << std::endl << std::endl;
    for(int i = 0; i < rotation_vector.size(); i++){
        std::cout << "Rotation:" << std::endl << rotation_vector[i] << std::endl;
        std::cout << "Translation:" << std::endl << translation_vector[i] << std::endl << std::endl;
    }

    /*--------------------------
    --------LAMBDA TWIST--------
    --------------------------*/

    cvl::Matrix<double, 3, 1> y0_cvl;
    cvl::Matrix<double, 3, 1> y1_cvl;
    cvl::Matrix<double, 3, 1> y2_cvl;
    y0_cvl = {triangle_in_world(0,0), triangle_in_world(1,0), triangle_in_world(2,0)};
    y1_cvl = {triangle_in_world(0,1), triangle_in_world(1,1), triangle_in_world(2,1)};
    y2_cvl = {triangle_in_world(0,2), triangle_in_world(1,2), triangle_in_world(2,2)};

    cvl::Matrix<double, 3, 1> x0_cvl;
    cvl::Matrix<double, 3, 1> x1_cvl;
    cvl::Matrix<double, 3, 1> x2_cvl;
    x0_cvl = {triangle_in_cam(0,0), triangle_in_cam(1,0), 1};
    x1_cvl = {triangle_in_cam(0,1), triangle_in_cam(1,1), 1};
    x2_cvl = {triangle_in_cam(0,2), triangle_in_cam(1,2), 1};

    cvl::Vector<cvl::Matrix<double,3,3>,4> Rs;
    cvl::Vector<cvl::Vector<double,3>,4> Ts;
    cvl::p3p_lambdatwist(y0_cvl, y1_cvl, y2_cvl, x0_cvl, x1_cvl, x2_cvl, Rs, Ts);

    std::cout << std::endl << "---LAMBDATWIST---" << std::endl;

    for(int i = 0; i<4; i++){
        std::cout << "Rotation:" << std::endl;
        std::cout << to_rodrigues(to_eigen_matrix(Rs(i))) << std::endl;
        
        std::cout << "Translation:" << std::endl;
        std::cout << Ts(i) << std::endl;
        std::cout << std::endl;
    }

}

int main(){
    
    
    test_p3p();


    return 0;
}
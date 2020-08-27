#include "Eigen/Dense"
#include "Eigen/Core"
#include <iostream>
#include <opencv2/opencv.hpp>

#define pi 3.14159265

template <class T>
Eigen::Matrix<T,3,3> make_rotation_matrix(const T& theta, const Eigen::Matrix<T,3,1>& w){
    Eigen::Matrix3d temp;
    temp = Eigen::AngleAxis<T>(theta, w);
    return temp;
}

int save_image(const Eigen::Matrix<double,3,3>& image_pts, const std::string& name){
    cv::Mat img(1280/2, 1920/2, CV_8UC3, cv::Scalar(255,255,255));
    if (img.empty()) 
    {
        std::cout << "\n Image not created. You"
                     " have done something wrong. \n";
        return -1;    // Unsuccessful.
    }

    for(int i=0; i<image_pts.cols(); i++){
        cv::circle(img,
                   cv::Point(image_pts(0,i)/2, image_pts(1,i)/2),
                   5,
                   cv::Scalar(0, 0, 255),
                   cv::FILLED);
    } 

    return 0;
}

bool is_in_frame(const Eigen::MatrixXd& pt,
                 const Eigen::Matrix<double,3,3>& k){
                     
    Eigen::Matrix<double,3,1> image_coords; 
    image_coords(0) = pt(0);
    image_coords(1) = pt(1);
    image_coords(2) = 1;
    image_coords = (k*image_coords).eval();

    if(image_coords(0) >= 1920 || image_coords(0) < 0){
        return false;
    }
    if(image_coords(1) >= 1080 || image_coords(1) < 0){
        return false;
    }
    return true;
}

template<class T, unsigned int rows, unsigned int cols>
Eigen::Matrix<T, rows, cols> gen_random_matrix(const T& min, const T& max){
    double range = max-min;

    Eigen::MatrixXd m = Eigen::MatrixXd::Random(rows,cols); 

    m = (m + Eigen::MatrixXd::Constant(rows,cols,1.))*range/2.; 
    m = (m + Eigen::MatrixXd::Constant(rows,cols,min));
    return m;
}

template<unsigned int cols>
Eigen::MatrixXd gen_pt_cloud(int num_pts=100, double midpoint=20, double max_dist=20){
    Eigen::MatrixXd pts = gen_random_matrix<double,4,cols>(midpoint-max_dist,midpoint+max_dist);
    pts.row(3) = Eigen::MatrixXd::Constant(1,pts.cols(),1.);
    return pts;
}

double fRand(const double& fMin, const double& fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

Eigen::Matrix<double, 4, 4> world_pts_in_cam(Eigen::MatrixXd& world_pts, 
                                             Eigen::MatrixXd& cam_pts, 
                                             const Eigen::Matrix<double,3,3>& k){
    bool run = true;

    Eigen::Vector3d w;
    Eigen::Matrix4d world_from_cam_transform;
    Eigen::MatrixXd cam_from_world_transform;
    Eigen::MatrixXd pts_in_cam_frame;
    Eigen::Vector4d world_from_cam_translation;
    
    w = gen_random_matrix<double, 3, 1> (-20, 20); // rotation axis
    w.normalize();
    
    world_from_cam_translation = gen_random_matrix<double, 4, 1> (-20, 0); //camera location
    world_from_cam_translation[3] = 1;

    Eigen::Vector3d cam_to_center = Eigen::MatrixXd::Constant(3,1,20) - world_from_cam_translation.block(0,0,3,1);
    cam_to_center.normalize();
    Eigen::Vector3d forward_vector;
    forward_vector << 0,0,1;
    w = forward_vector.cross(cam_to_center);
    w.normalize();

    double cos_angle = forward_vector.dot(cam_to_center);
    if(cos_angle < -1){
        cos_angle = -1;
    }
    else if(cos_angle > 1){
        cos_angle = 1;
    }

    world_from_cam_transform = Eigen::MatrixXd::Identity(4,4); 
    // world_from_cam_transform.block(0,0,3,3) = make_rotation_matrix(fRand(0, 2.0*pi), w);
    world_from_cam_transform.block(0,0,3,3) = make_rotation_matrix(acos(cos_angle), w);
    world_from_cam_transform.col(world_from_cam_transform.cols()-1) = world_from_cam_translation;
    cam_from_world_transform = world_from_cam_transform.inverse();

    pts_in_cam_frame = cam_from_world_transform*world_pts;

    cam_pts = pts_in_cam_frame.block(0,0,3,pts_in_cam_frame.cols());
    int i = 0;
    int num_cols = cam_pts.cols();
    while(i<num_cols){
        cam_pts(0,i) = cam_pts(0,i)/cam_pts(2,i);
        cam_pts(1,i) = cam_pts(1,i)/cam_pts(2,i);
        cam_pts(2,i) = 1;

        if(!is_in_frame(cam_pts.col(i), k)){
            removeColumn(cam_pts, i);
            removeColumn(world_pts, i);
            i--;
            num_cols = cam_pts.cols();
        }
        i++;
    }

    cam_pts = (k*cam_pts).eval();

    return cam_from_world_transform;
}


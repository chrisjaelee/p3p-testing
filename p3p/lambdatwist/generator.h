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

    cv::imwrite(name, img);

    return 0;
}

bool is_in_frame(const Eigen::MatrixXd& pts,
                 const Eigen::Matrix<double,3,3>& k,
                 const std::string& img_name = ""){
                     
    Eigen::Matrix<double,3,3> image_coords; 
    image_coords.block(0, 0, pts.rows(), pts.cols()) = pts;
    image_coords.row(2) << 1,1,1;
    image_coords = (k*image_coords).eval();

    for(int i=0; i<image_coords.cols(); i++){        
        if(image_coords(0,i) >= 1920 || image_coords(0,i) < 0){
            return false;
        }
        else if(image_coords(1,i) >= 1080 || image_coords(1,i) < 0){
            return false;
        }
    }

    if(img_name.size() > 0){
        save_image(image_coords, img_name);
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

double fRand(const double& fMin, const double& fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

Eigen::Matrix<double, 4, 4> world_pts_in_cam(Eigen::Matrix<double,4,3> world_pts, 
                      Eigen::MatrixXd &cam_pts, 
                      Eigen::Vector4d &world_from_cam_translation,
                      const Eigen::Matrix<double,3,3>& k){
    bool run = true;

    srand((unsigned int) time(0));

    Eigen::Vector3d w;
    Eigen::Matrix4d world_from_cam_transform;
    Eigen::Matrix4d cam_from_world_transform;
    Eigen::MatrixXd pts_in_cam_frame;
    
    while (run){
        w = gen_random_matrix<double, 3, 1> (-20, 20); // rotation axis
        w.normalize();
        
        world_from_cam_translation = gen_random_matrix<double, 4, 1> (-20, 20);
        world_from_cam_translation[3] = 1;

        world_from_cam_transform = Eigen::MatrixXd::Identity(4,4); 
        world_from_cam_transform.block(0,0,3,3) = make_rotation_matrix(fRand(0, 2.0*pi), w);
        world_from_cam_transform.col(world_from_cam_transform.cols()-1) = world_from_cam_translation;
        cam_from_world_transform = world_from_cam_transform.inverse();
        
        pts_in_cam_frame = cam_from_world_transform*world_pts;

        cam_pts.resize(2,3);
        
        cam_pts << pts_in_cam_frame(0,0), pts_in_cam_frame(0,1), pts_in_cam_frame(0,2),
                   pts_in_cam_frame(1,0), pts_in_cam_frame(1,1), pts_in_cam_frame(1,2);

        cam_pts.col(0) = cam_pts.col(0)/pts_in_cam_frame(2,0);
        cam_pts.col(1) = cam_pts.col(1)/pts_in_cam_frame(2,1);
        cam_pts.col(2) = cam_pts.col(2)/pts_in_cam_frame(2,2);

        // ADD CHECK IF IN FRAME HERE
        run = !is_in_frame(cam_pts, k);
        if (!run){
            return world_from_cam_transform.inverse();
        }
    }
}
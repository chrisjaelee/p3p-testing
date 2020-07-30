#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>

#define pi 3.14159265

template <class T>
void make_rotation_matrix(const T& theta, const Eigen::Matrix<T,3,1>& w, 
                          Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rotation){
    
    Eigen::Matrix3d temp;
    temp = Eigen::AngleAxis<T>(0, w);
    rotation = temp;
    
    rotation.conservativeResize(rotation.rows()+1, rotation.cols()+1);
    Eigen::Matrix<T,4,1> tmp;
    tmp << 0,0,0,1;

    rotation.row(rotation.rows()-1) = tmp;
    rotation.col(rotation.cols()-1) = tmp;
    

}

bool is_in_frame(const Eigen::MatrixXd& pts,
                 Eigen::Matrix<double,3,3> k){
                     
    Eigen::Matrix<double,3,3> image_coords; 
    image_coords.block(0, 0, pts.rows(), pts.cols()) = pts;
    image_coords.row(2) << 1,1,1;
    image_coords = (k*image_coords).eval();

    std::cout << std::endl;
    std::cout << "image_coordinates" << std::endl;
    std::cout << image_coords << std::endl;
    std::cout << std::endl;

    for(int i=0; i<image_coords.cols()-1; i++){
        if(image_coords(0,i) >= 1920 || image_coords(0,i) < 0){
            return false;
        }
        else if(image_coords(1,i) >= 1080 || image_coords(1,i) < 0){
            return false;
        }
    }
    return true;
}

void world_pts_in_cam(Eigen::Matrix<double,4,3> world_pts, 
                      Eigen::MatrixXd &cam_pts, 
                      Eigen::Vector4d world_from_cam_translation,
                      const Eigen::Matrix<double,3,3>& k,
                      const bool& verbose = false){
    bool run = true;
    // while (run){
        // creating rotation matrix
        Eigen::Vector3d w; // rotation axis
        w << 0, 0, 1;
        w.normalize();
        
        Eigen::MatrixXd world_from_cam_transform;
        make_rotation_matrix(pi, w, world_from_cam_transform);
        world_from_cam_transform.col(world_from_cam_transform.cols()-1) = world_from_cam_translation;
        Eigen::Matrix4d cam_from_world_transform = world_from_cam_transform.inverse();
        
        Eigen::MatrixXd pts_in_cam_frame;
        pts_in_cam_frame = cam_from_world_transform*world_pts;

        cam_pts.resize(2,3);
        
        cam_pts << pts_in_cam_frame(0,0), pts_in_cam_frame(0,1), pts_in_cam_frame(0,2),
                   pts_in_cam_frame(1,0), pts_in_cam_frame(1,1), pts_in_cam_frame(1,2);

        cam_pts.col(0) = cam_pts.col(0)/pts_in_cam_frame(2,0);
        cam_pts.col(1) = cam_pts.col(1)/pts_in_cam_frame(2,1);
        cam_pts.col(2) = cam_pts.col(2)/pts_in_cam_frame(2,2);

        // ADD CHECK IF IN FRAME HERE
        run = !is_in_frame(cam_pts, k);
        if (verbose && !run){
            std::cout << std::endl;
            std::cout << "Camera from world transform:" << std::endl;
            std::cout << cam_from_world_transform << std::endl;
            std::cout << std::endl;
        }
    // }
    
}
#include <iostream>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include "lambdatwist.p3p.h"
#include "utils/cvl/matrix.h"
#include <math.h>

template <class T, unsigned int rows, unsigned int cols>
Eigen::Matrix<T, rows, cols> cvlToEigen(const cvl::Matrix<T, rows, cols>& data){
    Eigen::Matrix<T, rows, cols> eigenMat;
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            eigenMat(i,j) = data(i,j);
        }
    }
    return eigenMat;
}

template <class T, unsigned int rows, unsigned int cols>
cvl::Matrix<T, rows, cols> eigenToCVL(const Eigen::Matrix<T, rows, cols>& data){
    cvl::Matrix<T, rows, cols> cvlMat;
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            cvlMat(i,j) = data(i,j);
        }
    }
    return cvlMat;
}

template <class T, unsigned int rows, unsigned int cols>
Eigen::Matrix<T, rows, cols> cvToEigen(const cv::Mat& data){
    Eigen::Matrix<T, rows, cols> eigenMatrix;
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            eigenMatrix(i,j) = data.at<T>(i,j);
        }
    }
    return eigenMatrix;
}

template<class T>
void rand3Pts(std::vector<cvl::Vector3<T>>& selectedWorldPts, std::vector<cvl::Vector3<T>>& selectedImagePts,
              const Eigen::Matrix<T, 3, Eigen::Dynamic>& worldPts, 
              const Eigen::Matrix<T, 3, Eigen::Dynamic>& imagePts){
    int matLen = worldPts.cols();

    std::vector<int> index;
    index.push_back(rand() % matLen);

    for(int i=0; i<3; i++){

        selectedWorldPts.push_back(eigenToCVL<T,3,1>(worldPts.col(index[i]) ));
        selectedImagePts.push_back(eigenToCVL<T,3,1>(imagePts.col(index[i]) ));
        
        index.push_back(rand() % matLen);
        for(int j=0; j<index.size()-1; j++){
            if(index[j] == index[index.size()-1]){
                index[index.size()-1] = rand() % matLen;
                j--;
            }
        }
    }
}

template<class T>
void reprojectPts(const Eigen::Matrix<T,3,Eigen::Dynamic>& worldPts, 
                        Eigen::Matrix<T,3,Eigen::Dynamic>& camPts, // output
                  const Eigen::Matrix<T,3,3>& rotation,
                  const Eigen::Matrix<T,3,1>& translation, // camera location
                  const Eigen::Matrix<T,3,3>& k){
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> ptsInCamFrame = worldPts;
    ptsInCamFrame.conservativeResize(4, worldPts.cols());
    for(int i=0;i<worldPts.cols();i++){
        ptsInCamFrame(3,i) = 1;
    }

    Eigen::Matrix<T,4,4> worldFromCamTransform; 
    worldFromCamTransform.block(0,0,3,3) = rotation;
    worldFromCamTransform.block(0,3,3,1) = translation;
    worldFromCamTransform.row(3) << 0,0,0,1;
    
    Eigen::Matrix<T,4,4> camFromWorldTransform = worldFromCamTransform.inverse();
    
    ptsInCamFrame = camFromWorldTransform*ptsInCamFrame;

    Eigen::MatrixXd temp_pt;
    camPts = ptsInCamFrame.block(0,0,3,ptsInCamFrame.cols());

    for(int i=0;i<camPts.cols();i++){
        camPts(0,i) = camPts(0,i)/ptsInCamFrame(2,i);
        camPts(1,i) = camPts(1,i)/ptsInCamFrame(2,i);
        camPts(2,i) = 1;
    }

    camPts = (k*camPts).eval();
}


template<class T>
T findError(const Eigen::Matrix<T,3,1>& estimate, const Eigen::Matrix<T,3,1>& actual){
    T distance = 0;
    for(int i=0;i<estimate.rows();i++){
        distance += pow(estimate(i)-actual(i),2);
    }
    return sqrt(distance);
}

template<class T>
int numInliers(const Eigen::Matrix<T,3,Eigen::Dynamic>& inputPts, 
               const Eigen::Matrix<T,3,Eigen::Dynamic>& reprojectedPts,
               const T& reprojectionError){
    int inliers = 0;

    for(int i=0;i<inputPts.cols();i++){

        
        T distance = 0;
        for(int j=0;j<2;j++){
            distance += pow(inputPts(j,i)-reprojectedPts(j,i),2);
        }
        distance = sqrt(distance);

        if(distance < reprojectionError){
            inliers++;
        }
    }

    return inliers;
}


template <class T>
Eigen::Matrix<T, 3, 1> to_rodrigues(const Eigen::Matrix<T, 3, 3>& data){
    Eigen::AngleAxis<double> aa (data);
    return aa.axis() * aa.angle();
}

template<class T>
void twistPnPRansac(const std::vector<cv::Point3_<T>>& worldPts, const std::vector<cv::Point_<T>>& imagePts, 
                    const cv::Mat& k, std::vector<cv::Mat>& r, std::vector<cv::Mat>& t,
                    int iterationsCount=100, T reprojectionError=8.0, const double& confidence=.99){
    srand(1234);

    Eigen::Matrix<T, 3, Eigen::Dynamic> eigenWorldPts; // fix initialization
    eigenWorldPts.resize(3,worldPts.size());
    Eigen::Matrix<T, 3, Eigen::Dynamic> eigenImagePts;
    eigenImagePts.resize(3,imagePts.size());
    Eigen::Matrix<T,3,3> eigenK = cvToEigen<T,3,3>(k);

    Eigen::Matrix<T,3,1> temp;
    for(int i=0; i<worldPts.size(); i++){
        temp << worldPts[i].x, 
                worldPts[i].y, 
                worldPts[i].z;
        eigenWorldPts.col(i) = temp;

        temp << imagePts[i].x, 
                imagePts[i].y,
                1;
        eigenImagePts.col(i) = eigenK.inverse()*temp;
    }

    std::vector<cvl::Vector3<T>> selectedWorldPts; 
    std::vector<cvl::Vector3<T>> selectedImagePts;

    cvl::Vector<cvl::Matrix<double,3,3>,4> Rs;
    cvl::Vector<cvl::Vector<double,3>,4> Ts;
    Eigen::Matrix<T,3,3> bestR;
    Eigen::Matrix<T,3,1> bestT;

    Eigen::Matrix<T,3,Eigen::Dynamic> reprojectedPts;
    Eigen::Matrix<T,3,3> rotation;
    Eigen::Matrix<T,3,1> translation;
    int bestInliers = 0;
    int inliers = 0;
    const int NUM_INNER_VALUES = 3;

    for(int iter=0; iter<iterationsCount; iter++){
        rand3Pts(selectedWorldPts, selectedImagePts, eigenWorldPts, eigenImagePts);
        
        // int valid = cvl::p3p_lambdatwist(selectedImagePts[0],selectedImagePts[1],selectedImagePts[2],
        //                                  selectedWorldPts[0],selectedWorldPts[1],selectedWorldPts[2],
        //                                  Rs, Ts);
        int valid = cvl::p3p_lambdatwist(eigenToCVL<T,3,1>(eigenImagePts.col(0)),
                                         eigenToCVL<T,3,1>(eigenImagePts.col(1)),
                                         eigenToCVL<T,3,1>(eigenImagePts.col(2)),
                                         eigenToCVL<T,3,1>(eigenWorldPts.col(0)),
                                         eigenToCVL<T,3,1>(eigenWorldPts.col(1)),
                                         eigenToCVL<T,3,1>(eigenWorldPts.col(2)),
                                         Rs, Ts);
        std::cout << "valid:\n" << valid << "\n\n";
        
        for(int i=0; i<valid; i++){
            rotation = cvlToEigen<T,3,3>(Rs[i]);
            std::cout << "rotation:\n" << to_rodrigues(rotation) << "\n\n";
            translation = cvlToEigen<T,3,1>(Ts[i]);
            reprojectPts(eigenWorldPts, reprojectedPts, 
                         rotation, translation, 
                         eigenK); 

            inliers = numInliers(eigenImagePts, reprojectedPts, reprojectionError);
            if(inliers > bestInliers){
                bestInliers = inliers;
                bestR = rotation;
                bestT = translation;
            }
        }
        
        T w = bestInliers/eigenImagePts.cols();
        T numIters = log10(1-confidence)/log10(1-pow(w,NUM_INNER_VALUES));
        if(numIters < iter){
            break;
        }

    }
    std::cout << bestR;
} 
cmake_minimum_required(VERSION 3.0)

find_package(OpenCV 3 REQUIRED)

find_package(Eigen3)
if(EIGEN3_FOUND)
    set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} ${VAR_SCOPE})
else()
    find_package(PkgConfig)
    pkg_check_modules(EIGEN3 eigen3)
    set(EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIRS} ${VAR_SCOPE})
endif()
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(p3p_tester main.cpp)
target_link_libraries(p3p_tester ${OpenCV_LIBRARIES})
target_include_directories(p3p_tester PUBLIC .)

PROJECT(test_main)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(SET CMP0015 NEW)

add_executable(IGCoupledMain IGCoupledMain.cpp)
target_link_libraries(IGCoupledMain IGCoupled)

find_package(Ceres REQUIRED)
target_link_libraries(IGCoupledMain ${CERES_LIBRARIES})

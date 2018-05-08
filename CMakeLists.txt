cmake_minimum_required(VERSION 2.8.3)

#project(manifoldReconstructor)  
project(filterByCamera)  
   
#----------------------- CMAKE MODULES ------------------------ 

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)  
  
# ----------------------- LIBRARIES ---------------------------
find_package(Boost COMPONENTS filesystem system serialization REQUIRED)
find_package(OpenCV REQUIRED)
find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)

# ----------------------- INCLUDES ----------------------------
include_directories(  ${EIGEN3_INCLUDE_DIR}
                      ${GMPXX_INCLUDE_DIR}
                      ${OpenCV_INCLUDE_DIRS} 
                      ${CGAL_INCLUDE_DIRS}
                      
                      ${PROJECT_SOURCE_DIR}/external/
                      ${PROJECT_SOURCE_DIR}/external/glm/glm/
                      ${PROJECT_SOURCE_DIR}/external/serialization/
                      ${PROJECT_SOURCE_DIR}/external/manifoldReconstructor/include/
                      ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/
                   
                      ${PROJECT_SOURCE_DIR}/include/coordinate_system_transform/
                      ${PROJECT_SOURCE_DIR}/include/utils/
                      
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/edge_managers/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/filtering/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/io/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/io/input/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/io/output/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/matching/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/matching/consensus_manager/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/matching/plg_matching/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/matching/polyline_matching/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/plgs/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/utils/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/utils/geometry/
                      ${PROJECT_SOURCE_DIR}/include/edgegraph3d/utils/globals/
                      ) 
# ----------------------- GCC FLAGS ----------------------------

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -w -O3 -lm -funroll-loops -std=c++11 -fpermissive -fopenmp -lboost_serialization")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -w -g -gdwarf-2 -g3 -lm -std=c++11 -fpermissive -fopenmp -lboost_serialization")

# ------------------------ MESSAGES ----------------------------

# ------------------------ BUILD -------------------------------

file(GLOB_RECURSE EdgeGraph3D_SOURCE 
                  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/*.cpp
                  ${PROJECT_SOURCE_DIR}/src/cam_parsers/*.cpp
		  		  ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/*.cpp
		  		  ${PROJECT_SOURCE_DIR}/external/manifoldReconstructor/src/*.cpp
                  )
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/driverForGraphClustering.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/convertFileToBinary.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/convertFileToEdgeList.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/driverForColoringExperiments.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/driverForRmat.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/driverForRGG.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/driverForPartitioningWithMetis.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/convertSnapFileToBinary.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/utilityGraphPartitioner.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/parallelLouvainMethodEager.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/driverProcessPowerGridNetworks.cpp)
list(REMOVE_ITEM EdgeGraph3D_SOURCE ${PROJECT_SOURCE_DIR}/external/grappolo-05-2014/driverForGraphClusteringBigMem.cpp)

file(GLOB_RECURSE json_to_ply_SOURCE 
				  ${PROJECT_SOURCE_DIR}/src/utils/json_to_ply.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/edge_graph_3d_utilities.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/io/output/output_sfm_data.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/edge_managers/segment_edge_manager.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/edge_managers/edge_manager.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/io/output/output_point_cloud.cpp
				  ${PROJECT_SOURCE_DIR}/external/manifoldReconstructor/src/OpenMvgParser.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/datatypes.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/geometry/geometric_utilities.cpp
                  )

file(GLOB_RECURSE coordinate_system_transform_SOURCE 
				  ${PROJECT_SOURCE_DIR}/src/utils/coordinate_system_transform.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/edge_graph_3d_utilities.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/io/output/output_sfm_data.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/edge_managers/segment_edge_manager.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/edge_managers/edge_manager.cpp
				  ${PROJECT_SOURCE_DIR}/external/manifoldReconstructor/src/OpenMvgParser.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/datatypes.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/geometry/geometric_utilities.cpp
				  ${PROJECT_SOURCE_DIR}/src/coordinate_system_transform/transform_coordinate_system.cpp
                  )
                  
file(GLOB_RECURSE filter_SOURCE 
				  ${PROJECT_SOURCE_DIR}/src/utils/filter.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/edge_graph_3d_utilities.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/io/output/output_sfm_data.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/edge_managers/segment_edge_manager.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/edge_managers/edge_manager.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/filtering/outliers_filtering.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/filtering/gauss_newton.cpp
				  ${PROJECT_SOURCE_DIR}/external/manifoldReconstructor/src/OpenMvgParser.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/datatypes.cpp
				  ${PROJECT_SOURCE_DIR}/src/edgegraph3d/utils/geometry/geometric_utilities.cpp
                  )
 
add_executable(EdgeGraph3D ${EdgeGraph3D_SOURCE})
add_executable(json_to_ply ${json_to_ply_SOURCE})
add_executable(coordinate_system_transform ${coordinate_system_transform_SOURCE})
add_executable(filter ${filter_SOURCE})

target_link_libraries(EdgeGraph3D 					${OpenCV_LIBS} 
                                            		${GMPXX_LIBRARIES}
                                            		${GMP_LIBRARIES}
                                            		${LAPACK_LIBRARIES}
                                            		${CGAL_LIBRARY}
                                            		${Boost_LIBRARIES}
                                            		${ASSIMP_LIBRARIES}  )
                                            
target_link_libraries(json_to_ply 					${OpenCV_LIBS} 
                                            		${GMPXX_LIBRARIES}
                                            		${GMP_LIBRARIES}
                                            		${Boost_LIBRARIES}
                                            		${CGAL_LIBRARY})     
                                            
target_link_libraries(coordinate_system_transform 	${OpenCV_LIBS} 
                                            		${GMPXX_LIBRARIES}
                                            		${GMP_LIBRARIES}
                                            		${Boost_LIBRARIES}
                                            		${CGAL_LIBRARY})      
                                            
target_link_libraries(filter 						${OpenCV_LIBS} 
                                            		${GMPXX_LIBRARIES}
                                            		${GMP_LIBRARIES}
                                            		${Boost_LIBRARIES}
                                            		${CGAL_LIBRARY})                                            

# ------------------------ INSTALLATION ------------------------


# --------------------------------------------------------------



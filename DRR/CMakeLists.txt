project(DRR)

find_package(ITK REQUIRED)
find_package(VTK REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${ITK_INCLUDE_DIRS})


add_executable(DRR DigitallyReconstructedRadiograph1.cxx)
target_link_libraries(DRR ${ITK_LIBRARIES} ${VTK_LIBRARIES})
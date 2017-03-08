option (CMAKE_BUILD_TYPE "CMake build type, e.g. RelWithDebInfo, Release, or Debug." RelWithDebInfo)

add_executable(
    pac_bam_scan
    pac_bam_scan.cpp)

#target_link_libraries (pac_bam_scan htslib)

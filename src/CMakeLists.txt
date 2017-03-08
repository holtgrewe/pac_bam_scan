option (CMAKE_BUILD_TYPE "CMake build type, e.g. RelWithDebInfo, Release, or Debug." RelWithDebInfo)

set (CMAKE_CXX_STANDARD 11)

set (ZLIB_INCLUDE_OPTS "")
set (ZLIB_LIBRARY_OPTS "")

find_package (ZLIB REQUIRED)
find_package (HTSlib REQUIRED)

include_directories (${HTSlib_INCLUDE_DIR})

add_executable(
    pac_bam_scan
    pac_bam_scan.cpp)

target_link_libraries (pac_bam_scan ${HTSlib_LIBRARY})

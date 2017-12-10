cmake_minimum_required(VERSION 3.8)
project(DOCIProject)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)
add_executable(DOCI_Head_run ${SOURCE_FILES})



include_directories(DOCI_Head)

add_subdirectory(DOCI_Head)
add_subdirectory(unit_tests)



target_link_libraries(DOCI_Head_run DOCI_Head)


find_package(eigen3)
if (EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIRS})
    target_link_libraries (DOCI_Head_run ${EIGEN3_LIBRARIES})
endif ()

find_package(Boost)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries (DOCI_Head_run ${Boost_LIBRARIES})
endif()

find_package(Armadillo)
if (Armadillo_FOUND)
    include_directories(${Armadillo_INCLUDE_DIR})
    target_link_libraries (DOCI_Head_run ${Armadillo_LIBRARIES})
endif()
find_package(hf 2.0.0 REQUIRED)
if (hf_FOUND)
    include_directories(${hf_INCLUDE_DIRS})
    target_link_libraries(DOCI_Head_run ${hf})
endif()

find_package(libwrp 2.1.1 REQUIRED)
if (libwrp_FOUND)
    include_directories(${libwrp_INCLUDE_DIRS})
    target_link_libraries(DOCI_Head_run ${libwrp_LIBRARIES})
endif()
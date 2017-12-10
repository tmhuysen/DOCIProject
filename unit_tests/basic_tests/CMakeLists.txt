include_directories(${googletest_SOURCE_DIR}/include ${googletest_SOURCE_DIR})

add_executable(basic_tests
        lol.cpp Doci_test.cpp Doci_real_test.cpp playtest.cpp)

target_link_libraries(basic_tests gtest gtest_main)
target_link_libraries(basic_tests DOCI_Head)

find_package(eigen3)
if (EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIRS})
    target_link_libraries (basic_tests ${EIGEN3_LIBRARIES})
endif ()

find_package(Boost)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries (basic_tests ${Boost_LIBRARIES})
endif()

find_package(Armadillo)
if (Armadillo_FOUND)
    include_directories(${Armadillo_INCLUDE_DIR})
    target_link_libraries (basic_tests ${Armadillo_LIBRARIES})
endif()
find_package(hf 2.0.0 REQUIRED)
if (hf_FOUND)
    include_directories(${hf_INCLUDE_DIRS})
    target_link_libraries(basic_tests ${hf_LIBRARIES})
endif()

find_package(libwrp 2.1.1 REQUIRED)
if (libwrp_FOUND)
    include_directories(${libwrp_INCLUDE_DIRS})
    target_link_libraries(basic_tests ${libwrp_LIBRARIES})
endif()
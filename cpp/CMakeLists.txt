cmake_minimum_required(VERSION 3.16)
project(qbp_capped_gbm LANGUAGES CXX)

add_library(qbp_capped_gbm src/capped_gbm.cpp)
target_compile_features(qbp_capped_gbm PUBLIC cxx_std_17)
target_include_directories(qbp_capped_gbm PUBLIC include)
target_compile_options(qbp_capped_gbm PRIVATE -Wall -Wextra -Wpedantic)

add_executable(qbp_gbm src/main.cpp)
target_link_libraries(qbp_gbm PRIVATE qbp_capped_gbm)

include(CTest)
if(BUILD_TESTING)
    add_executable(qbp_native_tests tests/native_tests.cpp)
    target_link_libraries(qbp_native_tests PRIVATE qbp_capped_gbm)
    add_test(NAME capped_gbm_native COMMAND qbp_native_tests)
endif()

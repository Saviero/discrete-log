cmake_minimum_required(VERSION 3.6)
project(discrete_log)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


set(SOURCE_FILES
        main.cpp
        FactorBaseGenerator.cpp
        FactorBaseGenerator.h
        PolynomialGenerator.cpp
        PolynomialGenerator.h
        SmoothGenerator.cpp
        SmoothGenerator.h
        )

add_executable(discrete_log ${SOURCE_FILES})

target_link_libraries(discrete_log ntl gmp armadillo)
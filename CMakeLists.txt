cmake_minimum_required(VERSION 3.0.0)

project(svd LANGUAGES CXX)

################################################################################

add_library(svd INTERFACE)

target_include_directories(svd INTERFACE src)

################################################################################

add_executable(test-svd EXCLUDE_FROM_ALL
	test/catch.cc
	test/calculate.cc
)

target_link_libraries(test-svd PRIVATE svd)

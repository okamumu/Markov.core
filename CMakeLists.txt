set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_minimum_required(VERSION 3.0.2)

enable_testing()

find_library(BLAS_LIBRARY blas)
SET(EXTRA_LIBS ${BLAS_LIBRARY})

# add_definitions("--std=c++11 -Wall")
add_definitions("-DF77BLAS -DF77LAPACK")
# add_definitions("-DDEBUG")

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wconversion -Wdisabled-optimization -Wendif-labels -Wfloat-equal -Winit-self -Winline -Wmissing-include-dirs -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wpacked -Wpointer-arith -Wredundant-decls -Wshadow -Wsign-promo -Wswitch-default -Wswitch-enum -Wvariadic-macros -Wwrite-strings")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -pg -DDEBUG")
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(test_dblas test/test_dblas.cpp)
target_link_libraries(test_dblas ${EXTRA_LIBS})
add_test(NAME dblas COMMAND $<TARGET_FILE:test_dblas>)
set_property(TEST dblas APPEND PROPERTY FAIL_REGULAR_EXPRESSION false)

add_executable(test_cppblas test/test_dblascpp.cpp)
target_link_libraries(test_cppblas ${EXTRA_LIBS})
add_test(NAME cppblas COMMAND $<TARGET_FILE:test_cppblas>)
set_property(TEST cppblas APPEND PROPERTY FAIL_REGULAR_EXPRESSION false)

add_executable(test_poisson test/test_poisson.cpp)
target_link_libraries(test_poisson ${EXTRA_LIBS})
add_test(NAME poisson COMMAND $<TARGET_FILE:test_poisson>)
set_property(TEST poisson APPEND PROPERTY FAIL_REGULAR_EXPRESSION false)

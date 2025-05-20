cmake_minimum_required(VERSION 3.16)
project(HodgkinHuxleyPlugin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Shared)

add_compile_options(-Wall -pedantic -Werror -g)
add_library(hodgkin_huxley SHARED "hodgkin_huxley.cpp")
target_include_directories(hodgkin_huxley PRIVATE ../src)
target_link_libraries(hodgkin_huxley PRIVATE ${CMAKE_DL_LIBS})

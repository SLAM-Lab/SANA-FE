cmake_minimum_required(VERSION 3.16)
project(HodgkinHuxleyPlugin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(CMAKE_BUILD_TYPE Shared)

add_compile_options(-Wall -pedantic -Werror -g)
add_library(hodgkin_huxley SHARED "hodgkin_huxley")
target_include_directories(hodgkin_huxley PRIVATE .)
#target_sources(hodgkinhuxley PRIVATE "hodgkin_huxley")
target_link_libraries(hodgkin_huxley PRIVATE ${CMAKE_DL_LIBS})

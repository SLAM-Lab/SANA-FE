cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(SANA-FE)

add_compile_options(
    -ldl -Wall -pedantic -Werror -g -fopenmp -Wno-stringop-truncation
)

file (GLOB SOURCE_FILES "*.cpp")
file (GLOB HEADER_FILES "*.hpp")

# Get the latest commit hash
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE)

add_compile_definitions(GIT_COMMIT="${GIT_COMMIT}")

find_package(pybind11 CONFIG)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

pybind11_add_module(simcpp 
	${SOURCE_FILES}
)
target_link_libraries(simcpp PRIVATE ${PYTHON_LIBRARIES})
target_link_libraries(simcpp PRIVATE pybind11::pybind11)
target_link_libraries(simcpp PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(simcpp PRIVATE ${CMAKE_DL_LIBS})
target_link_libraries(simcpp PUBLIC)

find_package(OpenMP REQUIRED)

add_executable(sim main.cpp)
target_sources(sim
    PRIVATE ${SOURCE_FILES}
)

target_link_libraries(sim PRIVATE ${PYTHON_LIBRARIES})
target_link_libraries(sim PRIVATE pybind11::pybind11)
target_link_libraries(sim PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(sim PRIVATE ${CMAKE_DL_LIBS})
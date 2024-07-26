cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
#set(PYTHON_BUILD_ENABLED True)
set(PYTHON_BUILD_ENABLED False)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(sana-fe)

add_compile_options(
    -Wall -pedantic -Werror -g -fopenmp
)

file (GLOB SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp")
file (GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")

# Get the latest commit hash
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE)

add_compile_definitions(GIT_COMMIT="${GIT_COMMIT}")

find_package(PythonLibs REQUIRED)
find_package(OpenMP REQUIRED)

#find_package(yaml-cpp 0.7.0 REQUIRED)
#### Try incorporating yaml-cpp using fetch content
include(FetchContent)

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 0.8.0
)
FetchContent_GetProperties(yaml-cpp)

if(NOT yaml-cpp_POPULATED)
  message(STATUS "Fetching yaml-cpp...")
  FetchContent_Populate(yaml-cpp)
  add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR})
endif()
############# End of yaml-cpp

############## rapid-yaml
#cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
#project(ryml-quickstart LANGUAGES CXX)
#
#if(RYML_BRANCH_NAME STREQUAL "")
#    find_program(GIT git REQUIRED)
#    execute_process(COMMAND git rev-parse --abbrev-ref HEAD
#        OUTPUT_VARIABLE RYML_BRANCH_NAME
#        OUTPUT_STRIP_TRAILING_WHITESPACE)
#endif()
#if(RYML_BRANCH_NAME STREQUAL HEAD)
#    set(RYML_BRANCH_NAME master)
#endif()
#message(STATUS "FetchContent from branch: ${RYML_BRANCH_NAME}")
#
#include(FetchContent)
#FetchContent_Declare(ryml
#    GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git
#    GIT_TAG ${RYML_BRANCH_NAME}
#    GIT_SHALLOW FALSE  # ensure submodules are checked out
#)
#FetchContent_MakeAvailable(ryml)

##add_executable(ryml-quickstart ../quickstart.cpp)
##target_link_libraries(ryml-quickstart ryml::ryml)

##add_custom_target(run ryml-quickstart
##    COMMAND $<TARGET_FILE:ryml-quickstart>
##    DEPENDS ryml-quickstart
##    COMMENT "running: $<TARGET_FILE:ryml-quickstart>")
##### END OF rapid-yaml

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${YAMLCPP_INCLUDE_DIRS})

# PyBind specific
if(PYTHON_BUILD_ENABLED)
find_package(pybind11 CONFIG)
pybind11_add_module(
    sanafecpp
	${SOURCE_FILES}
)
target_link_libraries(sanafecpp PRIVATE ${PYTHON_LIBRARIES})
target_link_libraries(sanafecpp PRIVATE pybind11::pybind11)
target_link_libraries(sanafecpp PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(sanafecpp PRIVATE ${CMAKE_DL_LIBS})
#target_link_libraries(sanafecpp PUBLIC yaml-cpp::yaml-cpp)
target_link_libraries(sanafecpp PRIVATE yaml-cpp)
endif()

#include(CMakePrintHelpers)
add_executable(sim main.cpp)
#cmake_print_variables(SOURCE_FILES)
list(FILTER SOURCE_FILES EXCLUDE REGEX "pymodule.cpp")
#cmake_print_variables(SOURCE_FILES)
target_sources(
    sim
    PRIVATE ${SOURCE_FILES}
)

#include(CMakePrintHelpers)
#cmake_print_variables(YAMLCPP_LIBRARIES)
if(PYTHON_BUILD_ENABLED)
target_link_libraries(sim PRIVATE ${PYTHON_LIBRARIES})
endif()
#end of PyBind specific
target_link_libraries(sim PRIVATE OpenMP::OpenMP_CXX)
#target_link_libraries(sim PRIVATE yaml-cpp)
target_link_libraries(sim PUBLIC yaml-cpp::yaml-cpp)
#target_link_libraries(sim PUBLIC ryml::ryml)
target_link_libraries(sim PRIVATE ${CMAKE_DL_LIBS})

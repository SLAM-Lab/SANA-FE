cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(POLICY CMP0148)
  cmake_policy(SET CMP0148 NEW)
endif()

# Set which targets to make, between a standalone executable and a Python lib
if (NOT DEFINED STANDALONE_BUILD_ENABLED)
set(STANDALONE_BUILD_ENABLED ON)
endif()

if (NOT DEFINED PYTHON_BUILD_ENABLED)
set(PYTHON_BUILD_ENABLED ON)
endif()

# Ensure that at least one build target is enabled
if(NOT STANDALONE_BUILD_ENABLED AND NOT PYTHON_BUILD_ENABLED)
    message(
        FATAL_ERROR
        "No build target enabled: Either STANDALONE_BUILD_ENABLED or PYTHON_BUILD_ENABLED must be ON."
    )
endif()

project(sana-fe)

add_compile_options(
    -Wall -pedantic -Werror -g -fopenmp -fPIC
)

include(CMakePrintHelpers)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
cmake_print_variables(SRC_DIR)
file (GLOB SOURCE_FILES "${SRC_DIR}/*.cpp")
file (GLOB HEADER_FILES "${SRC_DIR}/*.hpp")

cmake_print_variables(SOURCE_FILES)
cmake_print_variables(HEADER_FILES)

#file (GLOB SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp")
#file (GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")

# Get the latest commit hash
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE)

add_compile_definitions(GIT_COMMIT="${GIT_COMMIT}")

# Find Python using the new FindPython module
find_package(Python COMPONENTS Interpreter Development REQUIRED)
# Old version used PythonLibs
##find_package(PythonLibs REQUIRED)
find_package(OpenMP REQUIRED)

############## rapid-yaml
cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(ryml-quickstart LANGUAGES CXX)

if(RYML_BRANCH_NAME STREQUAL "")
    find_program(GIT git REQUIRED)
    execute_process(COMMAND git rev-parse --abbrev-ref HEAD
        OUTPUT_VARIABLE RYML_BRANCH_NAME
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
if(RYML_BRANCH_NAME STREQUAL HEAD)
    set(RYML_BRANCH_NAME master)
endif()
message(STATUS "FetchContent from branch: ${RYML_BRANCH_NAME}")

include(FetchContent)
FetchContent_Declare(ryml
    GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git
    GIT_TAG ${RYML_BRANCH_NAME}
    GIT_SHALLOW FALSE  # ensure submodules are checked out
)
FetchContent_MakeAvailable(ryml)
set(RYML_DEFAULT_CALLBACK_USES_EXCEPTIONS ON)
##### END OF rapid-yaml

include_directories(${PYTHON_INCLUDE_DIRS})

# PyBind specific
if(PYTHON_BUILD_ENABLED)
find_package(pybind11 CONFIG)
pybind11_add_module(
    sanafecpp
	${SOURCE_FILES}
)
target_link_libraries(sanafecpp PRIVATE ${PYTHON_LIBRARIES})
target_link_libraries(sanafecpp PRIVATE pybind11::pybind11)
target_link_libraries(sanafecpp PUBLIC ryml::ryml)
target_link_libraries(sanafecpp PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(sanafecpp PRIVATE ${CMAKE_DL_LIBS})
endif()

if(STANDALONE_BUILD_ENABLED)
add_executable(sim "${SRC_DIR}/main.cpp")
list(FILTER SOURCE_FILES EXCLUDE REGEX "pymodule.cpp")
target_sources(sim PRIVATE ${SOURCE_FILES})
endif()

if(STANDALONE_BUILD_ENABLED)
target_link_libraries(sim PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(sim PUBLIC ryml::ryml)
target_link_libraries(sim PRIVATE ${CMAKE_DL_LIBS})
endif()
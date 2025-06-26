cmake_minimum_required(VERSION 3.13)
include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
cmake_print_variables(CMAKE_BUILD_TYPE)

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
    -Wall -pedantic -Werror -g -fopenmp -fPIC -pthread
)

# Set default debug levels for tracing
if(NOT DEFINED DEBUG_LEVEL_ARCH)
    set(DEBUG_LEVEL_ARCH 0)
endif()

if(NOT DEFINED DEBUG_LEVEL_DESCRIPTION)
    set(DEBUG_LEVEL_DESCRIPTION 0)
endif()

if(NOT DEFINED DEBUG_LEVEL_MODELS)
    set(DEBUG_LEVEL_MODELS 0)
endif()

if(NOT DEFINED DEBUG_LEVEL_NET)
    set(DEBUG_LEVEL_NET 0)
endif()

if(NOT DEFINED DEBUG_LEVEL_PLUGINS)
    set(DEBUG_LEVEL_PLUGINS 0)
endif()

if(NOT DEFINED DEBUG_LEVEL_PYMODULE)
    set(DEBUG_LEVEL_PYMODULE 0)
endif()

if(NOT DEFINED DEBUG_LEVEL_CHIP)
    set(DEBUG_LEVEL_CHIP 0)
endif()

if(NOT DEFINED DEBUG_LEVEL_SCHEDULER)
    set(DEBUG_LEVEL_SCHEDULER 0)
endif()

# Validate debug levels
foreach(category ARCH CHIP DESCRIPTION MODELS NET PLUGINS PYMODULE SCHEDULER)
    if(NOT DEBUG_LEVEL_${category} MATCHES "^[0-3]$")
        message(FATAL_ERROR "DEBUG_LEVEL_${category} must be between 0 and 3")
    endif()
    # Add compile definitions for each category
    add_compile_definitions(DEBUG_LEVEL_${category}=${DEBUG_LEVEL_${category}})
endforeach()

# Set conditional defaults based on build type
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Release build defaults
    set(DEFAULT_STANDALONE_DEBUG ON)
    set(DEFAULT_STANDALONE_SOURCE OFF)  # No source info in release
    set(DEFAULT_PYTHON_DEBUG OFF)
    set(DEFAULT_PYTHON_SOURCE OFF)
else()
    # Debug/other build defaults
    set(DEFAULT_STANDALONE_DEBUG ON)
    set(DEFAULT_STANDALONE_SOURCE ON)
    set(DEFAULT_PYTHON_DEBUG ON)
    set(DEFAULT_PYTHON_SOURCE ON)
endif()

# Define options with conditional defaults
option(ENABLE_STANDALONE_DEBUG "Enable debug prints for standalone build" ${DEFAULT_STANDALONE_DEBUG})
option(ENABLE_SOURCE_INFO_STANDALONE "Include source info in standalone debug prints" ${DEFAULT_STANDALONE_SOURCE})
option(ENABLE_PYTHON_DEBUG "Enable debug prints for Python build" ${DEFAULT_PYTHON_DEBUG})
option(ENABLE_SOURCE_INFO_PYTHON "Include source info in standalone debug prints" ${DEFAULT_PYTHON_SOURCE})
cmake_print_variables(ENABLE_STANDALONE_DEBUG ENABLE_SOURCE_INFO_STANDALONE ENABLE_PYTHON_DEBUG ENABLE_SOURCE_INFO_PYTHON)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/plugins")
include(CMakePrintHelpers)
cmake_print_variables(SRC_DIR)
file (GLOB SOURCE_FILES "${SRC_DIR}/*.cpp")
file (GLOB HEADER_FILES "${SRC_DIR}/*.hpp")

cmake_print_variables(SOURCE_FILES)
cmake_print_variables(HEADER_FILES)

# Get the latest commit hash
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE)

add_compile_definitions(GIT_COMMIT="${GIT_COMMIT}")

# Find Python using the new FindPython module
option(PYTHON_FROM_SETUP "Building via setup.py" OFF)
set(PYTHON_VERSION "" CACHE STRING "Python version to find (e.g., 3.10, 3.11). Leave empty for system default.")
cmake_print_variables(PYTHON_VERSION)
if(PYTHON_BUILD_ENABLED)
    if(PYTHON_FROM_SETUP)
        if(NOT DEFINED PYTHON_EXECUTABLE OR NOT DEFINED PYTHON_INCLUDE_DIRS)
            message(FATAL_ERROR "PYTHON_EXECUTABLE and PYTHON_INCLUDE_DIRS must be provided when PYTHON_FROM_SETUP=ON")
        endif()
    else()
        if(PYTHON_VERSION)
            find_package(Python ${PYTHON_VERSION} EXACT COMPONENTS Interpreter Development REQUIRED)
        else()
            find_package(Python COMPONENTS Interpreter Development REQUIRED)
        endif()
    endif()
endif()

option(ENABLE_OPENMP "Enable OpenMP support for parallel processing" ON)
find_package(OpenMP)
find_package(Threads REQUIRED)

# OpenMP i.e. multithreaded builds are optional, but on by default
if (OpenMP_CXX_FOUND AND ENABLE_OPENMP)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION} and enabled")
else()
    if(NOT OpenMP_CXX_FOUND)
        message(STATUS "OpenMP not found - parallel code will be disabled")
    else()
        message(STATUS "OpenMP found but explicitly disabled - parallel code will be disabled")
    endif()
endif()

############## rapid-yaml
project(ryml-quickstart LANGUAGES CXX)
set(RYML_VERSION "v0.9.0")
message(STATUS "FetchContent for tag: ${RYML_VERSION}")

include(FetchContent)
FetchContent_Declare(ryml
    GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git
    GIT_TAG ${RYML_VERSION}
    GIT_SHALLOW FALSE  # ensure submodules are checked out
)
FetchContent_MakeAvailable(ryml)
set(RYML_DEFAULT_CALLBACK_USES_EXCEPTIONS ON)
##### END OF rapid-yaml

############## booksim2
project(booksim2-quickstart LANGUAGES CXX)
set(BOOKSIM_VERSION "110ad1b80e493241f6e57587bc11354ac84f91f8")
FetchContent_Declare(booksim
    GIT_REPOSITORY https://github.com/SLAM-Lab/booksim2-sanafe.git
    GIT_TAG ${BOOKSIM_VERSION}
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(booksim)
if(NOT booksim_POPULATED)
    FetchContent_Populate(booksim)
endif()

get_target_property(BOOKSIM_INCLUDE_DIRS booksim INTERFACE_INCLUDE_DIRECTORIES)
if(NOT BOOKSIM_INCLUDE_DIRS)
    # Fallback for older CMake versions
    set(BOOKSIM_INCLUDE_DIRS
        ${booksim_SOURCE_DIR}
        ${booksim_SOURCE_DIR}/src/arbiters
        ${booksim_SOURCE_DIR}/src/allocators
        ${booksim_SOURCE_DIR}/src/routers
        ${booksim_SOURCE_DIR}/src/networks
        ${booksim_SOURCE_DIR}/src/power
    )
    include_directories(${BOOKSIM_INCLUDE_DIRS})
endif()
#### END of booksim2

include_directories(${PYTHON_INCLUDE_DIRS})

cmake_print_variables(PyBind11_DIR)

# PyBind specific
if(PYTHON_BUILD_ENABLED)
    if(PYTHON_FROM_SETUP)
        # Force pybind11 to use the Python from setup.py
        set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
        set(Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
        set(PYBIND11_PYTHON_VERSION "") # Don't let pybind11 override
        message(STATUS "Forcing pybind11 to use Python: ${PYTHON_EXECUTABLE}")
    endif()

    find_package(pybind11 CONFIG)
    if (NOT PYTHON_EXECUTABLE)
        message("PYTHON_EXECUTABLE not set, default command: 'python3'")
        set(PYTHON_EXECUTABLE "python3")
    endif()
    if(NOT pybind11_FOUND OR pybind11_FOUND STREQUAL "0")
        message("PyBind11 not found, trying to set PYBIND_CMAKE_DIR using Python")
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
            OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        cmake_print_variables(PYBIND11_CMAKE_DIR)
        set(pybind11_DIR ${PYBIND11_CMAKE_DIR})
        find_package(pybind11 CONFIG REQUIRED)
    endif()
    pybind11_add_module(
        sanafecpp
        ${SOURCE_FILES}
    )

    target_link_libraries(sanafecpp PRIVATE ${PYTHON_LIBRARIES})
    target_link_libraries(sanafecpp PRIVATE pybind11::pybind11)
    target_link_libraries(sanafecpp PUBLIC ryml::ryml)
    target_link_libraries(sanafecpp PUBLIC booksim)
    target_link_libraries(sanafecpp PRIVATE ${CMAKE_DL_LIBS})
    if (OpenMP_CXX_FOUND AND ENABLE_OPENMP)
        target_link_libraries(sanafecpp PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(sanafecpp PRIVATE HAVE_OPENMP)
    endif()
    target_link_libraries(sanafecpp PRIVATE Threads::Threads)

    if(ENABLE_PYTHON_DEBUG)
        target_compile_definitions(sanafecpp PRIVATE ENABLE_DEBUG_PRINTS)
        if(ENABLE_SOURCE_INFO_PYTHON)
            target_compile_definitions(sanafecpp PRIVATE ENABLE_SOURCE_INFO)
        endif()
    endif()
endif()

if(STANDALONE_BUILD_ENABLED)
    add_executable(sim "${SRC_DIR}/main.cpp")
    list(FILTER SOURCE_FILES EXCLUDE REGEX "pymodule.cpp")
    target_sources(sim PRIVATE ${SOURCE_FILES})

    target_link_libraries(sim PUBLIC ryml::ryml)
    target_link_libraries(sim PUBLIC booksim)
    target_link_libraries(sim PRIVATE ${CMAKE_DL_LIBS})
    if (OpenMP_CXX_FOUND AND ENABLE_OPENMP)
        target_link_libraries(sim PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(sim PRIVATE HAVE_OPENMP)
    endif()
    target_link_libraries(sim PRIVATE Threads::Threads)

    if(ENABLE_STANDALONE_DEBUG)
        target_compile_definitions(sim PRIVATE ENABLE_DEBUG_PRINTS)
        if(ENABLE_SOURCE_INFO_STANDALONE)
            target_compile_definitions(sim PRIVATE ENABLE_SOURCE_INFO)
        endif()
    endif()
endif()

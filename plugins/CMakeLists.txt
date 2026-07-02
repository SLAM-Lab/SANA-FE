cmake_minimum_required(VERSION 3.16)
project(HodgkinHuxleyPlugin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Shared)

if(MSVC)
    add_compile_options(/W4 /WX /wd4127)  # Disable warn on empty control
else()
    add_compile_options(-Wall -pedantic -Werror)
endif()

add_library(hodgkin_huxley SHARED "hodgkin_huxley.cpp")
target_include_directories(hodgkin_huxley PRIVATE ../src)
target_link_libraries(hodgkin_huxley PRIVATE ${CMAKE_DL_LIBS})
set_target_properties(hodgkin_huxley PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(neurofem SHARED "neurofem.cpp")
target_include_directories(neurofem PRIVATE ../src)
target_link_libraries(neurofem PRIVATE ${CMAKE_DL_LIBS})
set_target_properties(neurofem PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

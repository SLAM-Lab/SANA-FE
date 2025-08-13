include(FetchContent)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# add unit_tests directory to the gtest suite
# MAKE SURE TO SET UP CMakeLists.txt IN ANY DIRECTORY YOU WANT TO ADD TO THE TEST SUITE
# you can follow the template in tests/unit_tests/CMakeLists.txt
add_subdirectory(unit_tests)

find_program(RUBY_EXE ruby REQUIRED)

# Same unique id + logs dir
execute_process(COMMAND git rev-parse --short HEAD
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
OUTPUT_VARIABLE GIT_SHA OUTPUT_STRIP_TRAILING_WHITESPACE)
string(TIMESTAMP TS "%Y%m%d-%H%M%S")
set(SANAFE_CI_ID "${TS}-${GIT_SHA}")

# Put logs in repo root
set(SANAFE_CI_LOG_DIR "${CMAKE_SOURCE_DIR}/logs/commit-${SANAFE_CI_ID}")
file(MAKE_DIRECTORY "${SANAFE_CI_LOG_DIR}")

# Register a Ruby checker as a CTest test
function(add_ruby_checker NAME SCRIPT_REL LABELS)
add_test(NAME ${NAME}
    COMMAND ${CMAKE_COMMAND} -E env
    SANAFE_CI_ID=${SANAFE_CI_ID}
    SANAFE_CI_LOG_DIR=${SANAFE_CI_LOG_DIR}
    ${RUBY_EXE} ${CMAKE_SOURCE_DIR}/${SCRIPT_REL}
)
set_tests_properties(${NAME} PROPERTIES
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    LABELS "${LABELS}"
)
endfunction()

add_ruby_checker(build    tests/ci/check_build.rb     "ci")
add_ruby_checker(format   tests/ci/check_format.rb    "dev;ci;lint")
add_ruby_checker(tidy     tests/ci/check_tidy.rb      "ci;static")
add_ruby_checker(cppcheck tests/ci/check_cppcheck.rb  "ci;static")
add_ruby_checker(dynamic  tests/ci/check_dynamic.rb   "ci;dynamic")
add_ruby_checker(perf     tests/ci/check_perf.rb      "ci;perf")

set_tests_properties(format   PROPERTIES DEPENDS build)
set_tests_properties(tidy     PROPERTIES DEPENDS format)
set_tests_properties(cppcheck PROPERTIES DEPENDS tidy)
set_tests_properties(dynamic  PROPERTIES DEPENDS cppcheck)
set_tests_properties(perf     PROPERTIES DEPENDS dynamic)

# Scenario targets: dev vs full
include(ProcessorCount)
ProcessorCount(NPROC)
if(NOT NPROC OR NPROC EQUAL 0)
set(NPROC 10)
endif()

add_custom_target(check-fast
DEPENDS unit_tests
COMMAND ${CMAKE_CTEST_COMMAND} -L dev --output-on-failure -j ${NPROC}
USES_TERMINAL)
add_custom_target(check-full
COMMAND ${CMAKE_CTEST_COMMAND} -L ci --output-on-failure -j ${NPROC}
USES_TERMINAL)
add_custom_target(check DEPENDS check-full)

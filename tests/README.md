# SANA-FE/tests
This directory contains all testing infrastructure and CI scripts for SANA-FE.

## Running Tests

### To run the full CI flow:
```bash
ruby tests/ci/run.rb
```

This will run the following scripts in ```tests/ci/```:
- ```check_build.rb```: verifies builds using GCC and Clang
- ```check_format.rb```: runs ```clang-format``` on all ```.cpp``` files in ```src/``` and ```plugins/```
- ```check_tidy.rb```: runs ```clang-tidy``` on all ```.cpp``` files in ```src/```
- ```check_cppcheck.rb```: runs ```cppcheck``` on all ```.cpp``` files in ```src/``` and ```plugins/```
- ```check_dynamic.rb```: builds project with testing enabled, then runs GoogleTest suite via CTest with Valgrind memory checks
- ```check_perf.rb```: runs an example scenario, records runtime, and issues a warning if runtime changed by more than 10%

### To use CTest directly, configure and build with CMake:
```bash
cmake -S . -B build -DENABLE_TESTING=ON
cmake --build build
```
Then, to run all CI-labeled tests:
```bash
cd build
ctest -L ci
```
Or, to run only format check and unit tests:
```bash
cd build
ctest -L dev
```
## Coverage Reports
To generate a code coverage report:
```bash
cd tests
./check_coverage.sh
```
This script will:
1. Clean and rebuild the project
2. Run all unit tests
3. Report coverage with gcovr

After it completes, a terminal summary will be printed and a full HTML report will be generated at ```tests/build/coverage.html```

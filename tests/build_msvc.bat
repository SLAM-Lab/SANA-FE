@echo off
setlocal

if "%VCVARS_PATH%"=="" (
  set "VCVARS_PATH=C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
)

REM %~dp0 is ...\tests\ ; the repo root is one level up.
pushd "%~dp0.."
if errorlevel 1 (
  echo Failed to enter script directory %~dp0
  exit /b 1
)

call "%VCVARS_PATH%"
if errorlevel 1 (
  echo Failed to load MSVC environment from "%VCVARS_PATH%"
  popd
  exit /b 1
)

cmake -S . -B build-msvc -G "Visual Studio 18 2026" -DCMAKE_BUILD_TYPE=Release -DPYTHON_BUILD_ENABLED=OFF -DENABLE_TESTING=ON
if errorlevel 1 ( popd & exit /b 1 )

cmake --build build-msvc --parallel 8
set BUILD_RC=%errorlevel%

popd
exit /b %BUILD_RC%

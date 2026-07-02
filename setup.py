"""Build script for sanafe's Python kernel."""
# pylint: disable=missing-class-docstring, missing-module-docstring
import os
import sys
import sysconfig
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(["cmake", "--version"])
        except OSError as exc:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions)
            ) from exc

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        print("Current directory:", os.getcwd())
        print("Source directory:", ext.sourcedir)
        print("External directory:", extdir)

        jobs = os.getenv('CMAKE_BUILD_PARALLEL_LEVEL', str(os.cpu_count() or 1))
        # Check for -j option
        if '-j' in sys.argv:
            # Find the index of '-j' and get the following number
            try:
                jobs_index = sys.argv.index('-j') + 1
                jobs = int(sys.argv[jobs_index])
                # Remove -j option and the following value from sys.argv
                sys.argv.pop(jobs_index)
                sys.argv.pop(jobs_index - 1)
            except (IndexError, ValueError):
                print("Warning: -j option requires a positive integer "
                      "argument, using default number of threads.")

        cmake_args = ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
                      "-DPYTHON_EXECUTABLE=" + sys.executable,
                      "-DPYTHON_INCLUDE_DIRS=" + sysconfig.get_path('include'),
                      "-DSTANDALONE_BUILD_ENABLED=OFF",
                      "-DPYTHON_BUILD_ENABLED=ON",
                      "-DPYTHON_FROM_SETUP=ON",
                      "-DBUILD_WHEEL=ON"]
        print(f"CMake Arguments: {cmake_args}")
        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg, "-j", jobs]

        if platform.system() == "Windows":
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
            cmake_args += [f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
            plat_name = self.plat_name
            if plat_name:
                plat_to_cmake = {
                    "win32": "Win32",
                    "win-amd64": "x64",
                    "win-arm64": "ARM64",
                }
                cmake_arch = plat_to_cmake.get(plat_name)
                if cmake_arch:
                    cmake_args += ["-A", cmake_arch]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]

        env = os.environ.copy()
        env["CXXFLAGS"] = (
            f"{env.get('CXXFLAGS', '')} "
            f"-DVERSION_INFO=\\'{self.distribution.get_version()}\\'"
        )
        env["CMAKE_BUILD_PARALLEL_LEVEL"] = str(jobs)
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        env["PYTHON_EXECUTABLE"] = sys.executable
        env["PYTHON_INCLUDE_DIRS"] = sysconfig.get_path('include')

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

setup(
    ext_modules=[CMakeExtension("sanafe")],
    cmdclass={"build_ext": CMakeBuild}
)

import os
import re
import sys
import sysconfig
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        print("Current directory:", os.getcwd())
        print("Source directory:", ext.sourcedir)
        print("External directory:", extdir)

        jobs = os.getenv('CMAKE_BUILD_PARALLEL_LEVEL', '1')  # Default to single-threaded build
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
                print("Warning: -j option requires a positive integer argument, using default number of threads.")

        cmake_args = ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
                      "-DPYTHON_EXECUTABLE=" + sys.executable,
                      "-DPYTHON_INCLUDE_DIRS=" + sysconfig.get_path('include'),
                      "-DSTANDALONE_BUILD_ENABLED=OFF",
                      "-DPYTHON_FROM_SETUP=ON"]
        print(f"CMake Arguments: {cmake_args}")
        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]

        env = os.environ.copy()
        env["CXXFLAGS"] = "{} -DVERSION_INFO=\\'{}\\'".format(env.get("CXXFLAGS", ""),
                                                              self.distribution.get_version())
        env["CMAKE_BUILD_PARALLEL_LEVEL"] = jobs
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

setup(
    name="sanafe",
    version="2.0.21",
    author="James Boyle",
    author_email="james.boyle@utexas.edu",
    description="SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/SLAM-Lab/SANA-FE",
    ext_modules=[CMakeExtension("sanafe")],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    python_requires=">=3.8",
    packages=find_packages()
)

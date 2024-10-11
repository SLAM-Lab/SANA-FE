import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion
import pybind11

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
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
        print(f"pybind directory: {pybind11.get_cmake_dir()}")
        cmake_args = ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
                      "-DPYTHON_EXECUTABLE=" + sys.executable,
                      f"-DPyBind11_DIR={pybind11.get_cmake_dir()}"]
        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = "{} -DVERSION_INFO=\\'{}\\'".format(env.get("CXXFLAGS", ""),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

setup(
    name="sanafe",
    version="0.0.1",
    author="James Boyle",
    author_email="james.boyle@utexas.edu",
    description="SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/SLAM-Lab/SANA-FE",
    ext_modules=[CMakeExtension("sanafe.cpp_module")],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=["pybind11>=2.6.0"],
    setup_requires=["pybind11>=2.6.0"],
    python_requires=">=3.6",
)
#!/usr/bin/env ruby

#TODO: Add support for MSVC builds

require 'fileutils'

commit_hash = `git rev-parse --short HEAD`.strip
log_dir = "logs/commit-#{commit_hash}"
FileUtils.mkdir_p(log_dir)

#gcc build
gcc_log = "#{log_dir}/build_gcc.log"
gcc_dir = "build"

puts "[GCC] Running CMake..."
gcc_cmake = system("cmake -S . -B #{gcc_dir} > #{gcc_log} 2>&1")

if gcc_cmake
  puts "[GCC] Running Make..."
  gcc_make = system("cmake --build #{gcc_dir} >> #{gcc_log} 2>&1")
end

if !gcc_cmake || !gcc_make
  puts "[GCC] Build FAILED. See #{gcc_log}"
  exit 1
else
  puts "[GCC] Build succeeded."
end

#clang 20 build
clang_bin = "/home/jab23579/useful/clang-20/bin/clang++"
clang_dir = "build_clang"
clang_log = "#{log_dir}/build_clang.log"
FileUtils.mkdir_p(clang_dir)

puts "[Clang] Running CMake with Clang 20..."
clang_cmake = system("cmake -S . -B #{clang_dir} -DCMAKE_CXX_COMPILER=#{clang_bin} > #{clang_log} 2>&1")

if clang_cmake
  puts "[Clang] Running Make with Clang 20..."
  clang_make = system("cmake --build #{clang_dir} >> #{clang_log} 2>&1")
end

if !clang_cmake || !clang_make
  puts "[Clang] Build FAILED. See #{clang_log}"
  exit 1
else
  puts "[Clang] Build succeeded."
end

exit 0



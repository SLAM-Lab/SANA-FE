#!/usr/bin/env ruby

require 'fileutils'

commit_hash = `git rev-parse --short HEAD`.strip
log_dir = "logs/commit-#{commit_hash}"
FileUtils.mkdir_p(log_dir)

#build standalone simulator
def build_cpp(label:, build_dir:, compiler: nil, log_file:)
  puts "[#{label}] Building standalone simulator..."

  cmake_cmd = "cmake -S . -B #{build_dir}"
  cmake_cmd += " -DCMAKE_CXX_COMPILER=#{compiler}" if compiler
  cmake_ok = system("#{cmake_cmd} > #{log_file} 2>&1")

  if cmake_ok
    make_ok = system("cmake --build #{build_dir} >> #{log_file} 2>&1")
  end

  if cmake_ok && make_ok
    puts "[#{label}] Simulator build: PASS"
    true
  else
    puts "[#{label}] Simulator build: FAIL (see #{log_file})"
    false
  end
end

#build python shared lib
def build_python(label:, build_dir:, compiler: nil, log_file:)
  puts "[#{label}] Building Python shared library..."

  cmake_cmd = "cmake -S . -B #{build_dir} -DENABLE_PYTHON=ON"
  cmake_cmd += " -DCMAKE_CXX_COMPILER=#{compiler}" if compiler
  cmake_ok = system("#{cmake_cmd} > #{log_file} 2>&1")

  if cmake_ok
    make_ok = system("cmake --build #{build_dir} >> #{log_file} 2>&1")
  end

  if cmake_ok && make_ok
    puts "[#{label}] Python build: PASS"
    true
  else
    puts "[#{label}] Python build: FAIL (see #{log_file})"
    false
  end
end

#clang path
clang_path = "/home/jab23579/useful/clang-20/bin/clang++"

#gcc builds
gcc_sim_ok = build_cpp(
  label: "GCC",
  build_dir: "build_gcc",
  log_file: "#{log_dir}/build_gcc_sim.log"
)
gcc_py_ok = build_python(
  label: "GCC",
  build_dir: "build_gcc_py",
  log_file: "#{log_dir}/build_gcc_py.log"
)

#clang builds
clang_sim_ok = build_cpp(
  label: "Clang",
  build_dir: "build_clang",
  compiler: clang_path,
  log_file: "#{log_dir}/build_clang_sim.log"
)
clang_py_ok = build_python(
  label: "Clang",
  build_dir: "build_clang_py",
  compiler: clang_path,
  log_file: "#{log_dir}/build_clang_py.log"
)

all_ok = gcc_sim_ok && gcc_py_ok && clang_sim_ok && clang_py_ok
puts all_ok ? "All builds succeeded." : "One or more builds failed."
exit(all_ok ? 0 : 1)



#!/usr/bin/env ruby

require 'fileutils'
require 'etc'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/cppcheck.log"
FileUtils.mkdir_p(log_dir)

venv_path   = File.expand_path("#{log_dir}/venv_cppcheck")
venv_python = File.join(venv_path, "bin", "python")

unless system("python3 -m venv #{venv_path}")
  puts "Failed to create virtualenv at #{venv_path}"
  exit(1)
end

setup_cmds = [
  "#{venv_python} -m pip install --upgrade pip wheel setuptools",
  "#{venv_python} -m pip install cmake ninja pybind11 scikit-build-core",
]
setup_cmds.each do |cmd|
  unless system("#{cmd} >> #{log_dir}/cppcheck_setup.log 2>&1")
    puts "ERROR: Failed to setup Python dependencies"
    exit(1)
  end
end

build_dir = "build_cppcheck"
FileUtils.rm_rf(build_dir)
FileUtils.mkdir_p(build_dir)
compile_commands = "#{build_dir}/compile_commands.json"

clang = ENV["CLANG_PATH"] || "clang++"
cmake_cmd = [
  "cmake", "-S .", "-B #{build_dir}",
  "-DCMAKE_CXX_COMPILER=#{clang}",
  "-DPYTHON_EXECUTABLE=#{venv_python}",
  "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
  "-DSTANDALONE_BUILD_ENABLED=OFF",
  "-DPYTHON_BUILD_ENABLED=ON",
  "-DENABLE_TESTING=OFF",
].join(" ")
unless system("#{cmake_cmd} > #{log_dir}/cppcheck_setup.log 2>&1")
  puts "ERROR: cmake configure failed"
  exit(1)
end

# Only check our own source — not _deps, not build artifacts
file_filters = %w[src/ plugins/].flat_map { |d| ["--file-filter=#{d}*"] }

cmd = [
  "cppcheck",
  "--project=#{compile_commands}",
  "--enable=warning,style,performance,portability",
  "--inline-suppr",
  "--suppress=missingIncludeSystem",
  "--suppress=*:*/_deps/*",
  "--suppress=*:*/build*/*",
  "--std=c++17",
  "--check-level=exhaustive",
  "-j#{Etc.nprocessors}",
  *file_filters,
  "--output-file=#{log_file}"
].join(" ")

puts "Running: #{cmd}"
success = system(cmd)

# cppcheck exits 0 even when it finds issues unless you ask it to fail.
# Treat any non-empty, non-trivial output as a failure.
issues = File.exist?(log_file) ? File.readlines(log_file).reject { |l| l.strip.empty? } : []

if issues.empty?
  puts "CPPCheck: PASS"
  exit(0)
else
  puts "CPPCheck: FAIL (#{issues.size} line(s), see #{log_file})"
  exit(1)
end

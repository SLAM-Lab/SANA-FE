#!/usr/bin/env ruby

require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/tidy.log"
FileUtils.mkdir_p(log_dir)

# Generate compile_commands.json if it doesn't exist
venv_path   = File.expand_path("#{File.dirname(log_file)}/venv_tidy")
venv_python = File.join(venv_path, "bin", "python")

unless system("python3 -m venv #{venv_path}")
  puts "[#{label}] Failed to create virtualenv at #{venv_path}"
  return false
end

# Quickly install the Python venv required for pybind11, to run CMake config
setup_cmds = [
    "#{venv_python} -m pip install --upgrade pip wheel setuptools",
    "#{venv_python} -m pip install cmake ninja pybind11 scikit-build-core"
  ]

  setup_cmds.each do |cmd|
    unless system("#{cmd} >> #{log_dir}/tidy_setup.log 2>&1")
      puts "ERROR: Failed to setup Python dependencies"
      return false
    end
  end

build_dir = "build_tidy"
compile_commands = "#{build_dir}/compile_commands.json"

unless File.exist?(compile_commands)
  puts "Generating compile_commands.json..."
  FileUtils.mkdir_p(build_dir)

  clang = ENV["CLANG_PATH"] || "clang++"
  cmake_cmd = [
    "cmake",
    "-S .",
    "-B #{build_dir}",
    "-DCMAKE_CXX_COMPILER=#{clang}",
    "-DPYTHON_EXECUTABLE=#{venv_python}",
    "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
  ].join(" ")

  unless system("#{cmake_cmd} > #{log_dir}/tidy_setup.log 2>&1")
    puts "ERROR: Failed to generate compile_commands.json"
    exit(1)
  end
end

cpp_files = Dir.glob("src/*.cpp")
failed_files = []

File.open(log_file, "w") do |log|
  cpp_files.each do |file|
    log.puts "----- #{file} -----"
    result = `clang-tidy -p #{build_dir} #{file} -- -I./src 2>&1`
    log.puts result

    if result.include?("warning:") || result.include?("error:")
      failed_files << file
    end
  end

  if failed_files.empty?
    log.puts "\nAll files passed clang-tidy check."
  else
    log.puts "\n#{failed_files.size} file(s) had warnings/errors."
  end
end

if failed_files.empty?
  puts "Tidy Check: PASS"
else
  puts "Tidy Check: FAIL (#{failed_files.size} file(s) had issues, see #{log_file})"
end

exit(failed_files.empty? ? 0 : 1)

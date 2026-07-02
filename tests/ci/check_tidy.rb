#!/usr/bin/env ruby

require 'fileutils'
require 'etc'
require 'open3'
require 'thread'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/tidy.log"
FileUtils.mkdir_p(log_dir)

venv_path   = File.expand_path("#{log_dir}/venv_tidy")
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
  unless system("#{cmd} >> #{log_dir}/tidy_setup.log 2>&1")
    puts "ERROR: Failed to setup Python dependencies"
    exit(1)
  end
end

build_dir = "build_tidy"
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

unless system("#{cmake_cmd} > #{log_dir}/tidy_setup.log 2>&1")
  puts "ERROR: cmake configure failed"
  exit(1)
end

# Decide what to check
target = ENV["TIDY_TARGET"]  # e.g. "src/main.cpp" to check just one file
cpp_files = if target
              [target]
            else
              Dir.glob("src/*.cpp")
            end

jobs = (ENV["TIDY_JOBS"] || Etc.nprocessors).to_i
queue = Queue.new
cpp_files.each { |f| queue << f }

results = {}        # file => { output:, status: }
results_mutex = Mutex.new
diag_re = /^[^:\n]+:\d+:\d+:\s+(warning|error):/

total = cpp_files.size
done = 0
progress_mutex = Mutex.new

workers = jobs.times.map do
  Thread.new do
    while (file = queue.pop(true) rescue nil)
      # Open3.capture2e merges stdout+stderr into one string, no shell needed
      output, status = Open3.capture2e("clang-tidy", "-p", build_dir,
                                       "--quiet", file)
      results_mutex.synchronize { results[file] = { output: output, status: status } }
      progress_mutex.synchronize do
        done += 1
        $stderr.print "\r[#{done}/#{total}] checked"
      end
    end
  end
end

workers.each(&:join)
$stderr.puts

# Write results in deterministic (input) order
failed_files = []
File.open(log_file, "w") do |log|
  cpp_files.each do |file|
    r = results[file]
    log.puts "----- #{file} -----"
    log.puts r[:output]
    failed_files << file if r[:output] =~ diag_re
  end

  log.puts
  if failed_files.empty?
    log.puts "All files passed clang-tidy check."
  else
    log.puts "#{failed_files.size} file(s) had warnings/errors:"
    failed_files.each { |f| log.puts "  #{f}" }
  end
end

if failed_files.empty?
  puts "Tidy Check: PASS"
  exit(0)
else
  puts "Tidy Check: FAIL (#{failed_files.size} file(s), see #{log_file})"
  exit(1)
end

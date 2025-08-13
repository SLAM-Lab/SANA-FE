#!/usr/bin/env ruby

require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
build_log_file = "#{log_dir}/dynamic_build.log"
test_log_file = "#{log_dir}/dynamic_test.log"
FileUtils.mkdir_p(log_dir)

# timeout after 30 seconds
timeout = 30

label = "gtest"

puts "\n[#{label}] Running dynamic tests..."

puts "[#{label}] Building the project with testing enabled..."

cmake = system("cmake -DENABLE_TESTING=ON -DPYTHON_BUILD_ENABLED=OFF -DSTANDALONE_BUILD_ENABLED=OFF -DCMAKE_BUILD_TYPE=Debug -S . -B build > #{build_log_file} 2>&1")

if !cmake
  puts "[#{label}] CMake configuration failed. See #{build_log_file} for details."
  puts "[#{label}] Dynamic Tests: FAIL"
  exit 3
end

build = system("cmake --build build -j 10 >> #{build_log_file} 2>&1")

if !build
  puts "[#{label}] Build failed. See #{build_log_file} for details."
  puts "[#{label}] Dynamic Tests: FAIL"
  exit 3
end

puts "[#{label}] Build successful. Running tests..."

test = system("ctest --memcheck --test-dir build --output-on-failure --timeout #{timeout} > #{test_log_file} 2>&1")
exitcode = $?.exitstatus

if exitcode != 0
  puts "[#{label}] One or more tests failed. See #{test_log_file} for details."
  puts "[#{label}] Dynamic Tests: FAIL"
  exit 1
end

puts "[#{label}] Tests completed successfully, now checking for memory leaks..."

build_dir = "build"
mem_log_dir = "#{build_dir}/Testing/Temporary"
pattern = "MemoryChecker.*\\.log$"
summary_log = File.join(mem_log_dir, "LastMemCheck.log")

unless Dir.exist?(mem_log_dir)
  puts "[#{label}] Memory log directory #{mem_log_dir} does not exist. Something went wrong."
  exit 4
end

log_files = Dir.entries(mem_log_dir).select { |f| f.match?(pattern) }
puts "[#{label}] Found #{log_files.size} Valgrind log(s)."

leaks_found = false

log_files.each do |filename|
  path = File.join(mem_log_dir, filename)
  contents = File.read(path)

  if contents.include?("definitely lost")
    puts "[#{label}] Leak detected in #{filename}"
    leaks_found = true
  end
end

if leaks_found
  puts "[#{label}] Memory leaks detected in one or more Valgrind logs. See #{mem_log_dir} for details."
  puts "[#{label}] Dynamic Tests: FAIL"
  exit 2
else
  puts "[#{label}] All tests passed successfully with no memory leaks. Hooray!"
  puts "[#{label}] Dynamic Tests: PASS"
  exit 0
end
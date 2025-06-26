#!/usr/bin/env ruby

require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
build_log_file = "#{log_dir}/dynamic_build.log"
test_log_file = "#{log_dir}/dynamic_test.log"
FileUtils.mkdir_p(log_dir)

# timeout after 30 seconds
timeout = 30

puts "Running dynamic tests..."

puts "Building the project with testing enabled..."

cmake = system("cmake -DENABLE_TESTING=ON -DPYTHON_BUILD_ENABLED=OFF -DSTANDALONE_BUILD_ENABLED=OFF -DCMAKE_BUILD_TYPE=Debug -S . -B build > #{build_log_file} 2>&1")

if !cmake
  puts "CMake configuration failed. See #{build_log_file} for details."
  puts "Dynamic Tests: FAIL"
  exit 3
end

build = system("cmake --build build -j 10 >> #{build_log_file} 2>&1")

if !build
  puts "Build failed. See #{build_log_file} for details."
  puts "Dynamic Tests: FAIL"
  exit 3
end

puts "Build successful. Running tests..."

test = system("ctest --memcheck --test-dir build --output-on-failure --timeout #{timeout} > #{test_log_file} 2>&1")
exitcode = $?.exitstatus

puts "Tests completed, now checking for memory leaks..."

build_dir = "build"
mem_log_dir = "#{build_dir}/Testing/Temporary"
pattern = "MemoryChecker.*\\.log$"
summary_log = File.join(mem_log_dir, "LastMemCheck.log")

unless Dir.exist?(mem_log_dir)
  puts "Memory log directory #{mem_log_dir} does not exist. Something went wrong."
  exit 2
end

log_files = Dir.entries(mem_log_dir).select { |f| f.match?(pattern) }
puts "Found #{log_files.size} Valgrind log(s)."

leaks_found = false

log_files.each do |filename|
  path = File.join(mem_log_dir, filename)
  contents = File.read(path)

  if contents.include?("definitely lost")
    puts "Leak detected in #{filename}"
    leaks_found = true
  end
end

if leaks_found
  puts "Memory leaks detected in one or more Valgrind logs. See #{mem_log_dir} for details."
  puts "Dynamic Tests: FAIL"
  exit 2
end
if exitcode != 0
    puts "One or more tests failed. See #{test_log_file} for details."
    puts "Dynamic Tests: FAIL"
    exit 1
else
    puts "All tests passed successfully. Hooray! Check build/Testing/Temporary for memory leak reports."
    puts "Dynamic Tests: PASS"
    exit 0
end
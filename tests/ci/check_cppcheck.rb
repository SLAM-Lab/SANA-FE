#!/usr/bin/env ruby

require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/cppcheck.log"
FileUtils.mkdir_p(log_dir)

cpp_files = Dir.glob("src/*.cpp") + Dir.glob("plugins/*.cpp")
failed = false

File.open(log_file, "w") do |log|
  log.puts "Running CPPCheck..."
  cpp_files.each do |file|
    log.puts "----- #{file} -----"
    result = `cppcheck --enable=all --inconclusive --quiet --std=c++20 --language=c++ --suppress=missingIncludeSystem #{file} 2>&1`
    log.puts result
    failed ||= !result.strip.empty?
  end
end

puts failed ? "CPPCheck: FAIL (see #{log_file})" : "CPPCheck: PASS"
exit(failed ? 1 : 0)

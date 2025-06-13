#!/usr/bin/env ruby

require 'fileutils'

#paths
commit_hash = `git rev-parse --short HEAD`.strip
log_file = "logs/commit-#{commit_hash}/format.log"

#find cpp files
cpp_files = Dir.glob("**/*.cpp")

failed_files = []

File.open(log_file, "w") do |log|
  cpp_files.each do |file|
    result = `clang-format --dry-run --Werror #{file} 2>&1`
    if $?.exitstatus != 0
      failed_files << file
      log.puts "[FAIL] #{file}"
      log.puts result
    else
      log.puts "[PASS] #{file}"
    end
  end

  if failed_files.empty?
    log.puts "\nAll files passed clang-format check."
  else
    log.puts "\n#{failed_files.size} file(s) failed formatting."
  end
end

exit(failed_files.empty? ? 0 : 1)

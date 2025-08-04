#!/usr/bin/env ruby
#TODO: break format
require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/format.log"
FileUtils.mkdir_p(log_dir)

#find cpp files
cpp_files = Dir.glob("src/*.cpp") + Dir.glob("plugins/*.cpp")

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

if failed_files.empty?
  puts "Format Check: PASS"
else
  puts "Format Check: FAIL (#{failed_files.size} file(s) failed, see #{log_file})"
end

exit(failed_files.empty? ? 0 : 1)

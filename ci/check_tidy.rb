#!/usr/bin/env ruby

require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/tidy.log"
FileUtils.mkdir_p(log_dir)

cpp_files = Dir.glob("src/*.cpp")  #
failed_files = []

File.open(log_file, "w") do |log|
  cpp_files.each do |file|
    log.puts "----- #{file} -----"
    result = `clang-tidy #{file} -- -I./src 2>&1`
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

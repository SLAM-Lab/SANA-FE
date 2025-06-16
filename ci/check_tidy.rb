#!/usr/bin/env ruby

require 'fileutils'

commit_hash = `git rev-parse --short HEAD`.strip
log_file = "logs/commit-#{commit_hash}/tidy.log"
FileUtils.mkdir_p("logs/commit-#{commit_hash}")

#all .cpp files in src
cpp_files = Dir.glob("src/**/*.cpp")

failed = false

File.open(log_file, "w") do |log|
  cpp_files.each do |file|
    log.puts "----- #{file} -----"
    result = `clang-tidy #{file} -- -I./src 2>&1`
    log.puts result

    failed ||= result.include?("warning:") || result.include?("error:")
  end
end

exit failed ? 1 : 0

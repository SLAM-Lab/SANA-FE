#!/usr/bin/env ruby

require 'fileutils'

#paths
commit_hash = `git rev-parse --short HEAD`.strip
log_dir = "logs/commit-#{commit_hash}"
FileUtils.mkdir_p(log_dir)

puts "Running SANA-FE CI for commit #{commit_hash}"
puts "Logs will be saved to #{log_dir}/"
puts "-------------------------------"

#run build checker
build_status = system("ruby ci/check_build.rb")

#run clang-format
format_status = system("ruby ci/check_format.rb")

#run clang-tidy
tidy_status = system("ruby ci/check_tidy.rb")

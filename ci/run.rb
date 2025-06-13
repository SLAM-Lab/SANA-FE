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
puts "Build Check: #{build_status ? 'PASS' : 'FAIL'}"

#run format checker
format_status = system("ruby ci/check_format.rb")
puts "Format Check: #{format_status ? 'PASS' : 'FAIL'}"
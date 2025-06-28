#!/usr/bin/env ruby

require 'fileutils'

timestamp = Time.now.strftime("%Y%m%d-%H%M%S")            
commit_hash = `git rev-parse --short HEAD`.strip          
unique_id = "#{timestamp}-#{commit_hash}"                 
log_dir = "logs/commit-#{unique_id}"   
FileUtils.mkdir_p(log_dir)

ENV["SANAFE_CI_LOG_DIR"] = log_dir
ENV["SANAFE_CI_ID"] = unique_id

puts "Running SANA-FE CI for commit #{commit_hash}"
puts "Logs will be saved to #{log_dir}/"
puts "-------------------------------"

build_status = system("ruby ci/check_build.rb")    
format_status = system("ruby ci/check_format.rb")   
tidy_status = system("ruby ci/check_tidy.rb")       
cppcheck_status = system("ruby ci/check_cppcheck.rb")
#TODO: Optional coverity integration
dynamic_status = system("ruby ci/check_dynamic.rb")
perf_status = system("ruby ci/check_perf.rb")

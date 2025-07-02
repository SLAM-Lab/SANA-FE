#!/usr/bin/env ruby

require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/perf.log"
FileUtils.mkdir_p(log_dir)

#test parameters
sim_path = "./sim"
arch_file = "arch/example.yaml"
snn_file = "snn/example.yaml"
steps = 100000

perf_output = `(/usr/bin/time -f "%e" #{sim_path} #{arch_file} #{snn_file} #{steps} > /dev/null) 2>&1`
runtime = perf_output.to_f

baseline_file = "ci/perf_baseline.txt"
commit_perf_file = "#{log_dir}/perf_runtime.txt"

#first time 
if !File.exist?(baseline_file)
  File.write(baseline_file, "#{runtime.round(3)}\n")
  File.write(log_file, <<~LOG)
    Baseline not found. Created a new baseline:
    Baseline Runtime: #{runtime.round(3)} sec

    To reset baseline manually later:
      cp #{commit_perf_file} #{baseline_file}
  LOG
  puts "Performance Check: Baseline Created (#{runtime.round(3)} sec)"
  File.write(commit_perf_file, "#{runtime.round(3)}\n")
  exit 0
end

#compare with baseline
baseline = File.read(baseline_file).to_f
delta = (runtime - baseline).abs / baseline
warning = delta > 0.10

#write runtime
File.write(commit_perf_file, "#{runtime.round(3)}\n")

#write to log
File.open(log_file, "w") do |f|
  f.puts "Baseline Runtime: #{baseline.round(3)} sec"
  f.puts "Current Runtime:  #{runtime.round(3)} sec"
  f.puts "Delta:           #{(delta * 100).round(2)}%"
  f.puts warning ? "Performance Check: WARNING (runtime changed > 10%)" : "Performance Check: PASS"
end

#results
puts warning ? "Performance Check: WARNING (see #{log_file})" : "Performance Check: PASS"
exit(warning ? 1 : 0)

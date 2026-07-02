#!/usr/bin/env ruby

require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
FileUtils.mkdir_p(log_dir)

# Path to vcvars64.bat as a *Linux* env var. We translate it to a Windows
# path with wslpath before handing it to cmd.exe. Falls back to the common
# VS 2026 Community location if unset.
DEFAULT_VCVARS =
  "C:\\Program Files\\Microsoft Visual Studio\\18\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat"

def to_windows_path(p)
  # If it already looks like a Windows path (has a drive letter or backslashes),
  # leave it alone; otherwise translate from a /mnt/c-style Linux path.
  return p if p =~ /\A[A-Za-z]:\\/ || p.include?("\\")
  out = `wslpath -w #{p.shellescape} 2>/dev/null`.strip
  out.empty? ? p : out
end

require 'shellwords'

def build_msvc(label:, log_file:)
  puts "[#{label}] Building standalone simulator (MSVC)..."
  FileUtils.mkdir_p(File.dirname(log_file))

  vcvars = ENV["VCVARS_PATH"]
  vcvars = to_windows_path(vcvars) if vcvars
  vcvars ||= DEFAULT_VCVARS

  # Hand VCVARS_PATH to the batch script through cmd.exe's environment.
  # cmd.exe /c runs build_msvc.bat on the Windows host; WSL interop bridges it.
  bat_path = File.expand_path("../build_msvc.bat", __dir__)
  bat_win  = to_windows_path(bat_path)   # -> \\wsl.localhost\Ubuntu\...\build_msvc.bat
  cmd_exe  = ENV["CMD_EXE"] || "/mnt/c/Windows/System32/cmd.exe"
  cmd = %(VCVARS_PATH=#{vcvars.shellescape} #{cmd_exe.shellescape} /c #{bat_win.shellescape})

  ok = false
  File.open(File.expand_path(log_file), "w") do |log|
    log.sync = true
    log.puts("$ #{cmd}")
    IO.popen(cmd + " 2>&1") { |io| io.each { |line| log.write(line) } }
    ok = $?.success?
  end

  if ok
    puts "[#{label}] MSVC build: PASS"
  else
    puts "[#{label}] MSVC build: FAIL (see #{log_file})"
  end
  ok
end

msvc_ok = build_msvc(
  label: "MSVC",
  log_file: "#{log_dir}/build_msvc_sim.log"
)

puts msvc_ok ? "MSVC build succeeded." : "MSVC build failed."
exit(msvc_ok ? 0 : 1)

#!/usr/bin/env ruby

require 'fileutils'

log_dir = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/build.log"
FileUtils.mkdir_p(log_dir)


#build standalone sim
def build_cpp(label:, build_dir:, compiler: nil, log_file:)
  puts "[#{label}] Building standalone simulator..."

  cmake_cmd = "cmake -S . -B #{build_dir}"  #construct cmake config command using source and build dirs
  cmake_cmd += " -DCMAKE_CXX_COMPILER=#{compiler}" if compiler  #add compiler option if provided
  cmake_ok = system("#{cmake_cmd} > #{log_file} 2>&1")  #run cmake, direct output to log file

  if cmake_ok
    make_ok = system("cmake --build #{build_dir} --parallel 8 >> #{log_file} 2>&1")  #if cmake works, build and append output to log
  end

  if cmake_ok && make_ok
    puts "[#{label}] Simulator build: PASS"
    true
  else
    puts "[#{label}] Simulator build: FAIL (see #{log_file})"
    false
  end
end

#build and install python .so
def build_python(label:, build_dir:, log_file:)
  puts "[#{label}] Building python extension..."

  FileUtils.rm_f("CMakeCache.txt")
  FileUtils.mkdir_p(build_dir)
  FileUtils.mkdir_p(File.dirname(log_file))

  full_log_path = File.expand_path(log_file)

  #create venv name
  unique_id = ENV["SANAFE_CI_ID"]
  venv_path = File.expand_path("#{File.dirname(log_file)}/venv")
  venv_python = File.join(venv_path, "bin", "python")

  FileUtils.mkdir_p("venvs")
  unless system("python3 -m venv #{venv_path}")
    puts "[#{label}] Failed to create virtualenv at #{venv_path}"
    return false
  end

  install_ok = false
  import_ok = false

  File.open(full_log_path, "w") do |log|
    Dir.chdir(build_dir) do
      IO.popen("#{venv_python} -m pip install .. > /dev/null 2>&1") do |io|
        io.each { |line| log.puts line }
      end
      install_ok = $?.success?
    end
  end

puts "[#{label}] Python build: #{install_ok ? 'PASS' : "FAIL (see #{log_file})"}"
  install_ok
end

#set clang and gcc path
clang_path = ENV["CLANG_PATH"]
gcc_path = ENV["GCC_PATH"]

#gcc builds
gcc_sim_ok = build_cpp(
  label: "GCC",
  build_dir: "build_gcc",
  compiler: gcc_path,
  log_file: "#{log_dir}/build_gcc_sim.log"
)
gcc_py_ok = build_python(
  label: "GCC",
  build_dir: "build_gcc_py",
  log_file: "#{log_dir}/build_gcc_py.log"
)

#clang builds
clang_sim_ok = build_cpp(
  label: "Clang",
  build_dir: "build_clang",
  compiler: clang_path,
  log_file: "#{log_dir}/build_clang_sim.log"
)
clang_py_ok = build_python(
  label: "Clang",
  build_dir: "build_clang_py",
  log_file: "#{log_dir}/build_clang_py.log"
)

#final summary
all_ok = gcc_sim_ok && gcc_py_ok && clang_sim_ok && clang_py_ok
puts all_ok ? "All builds succeeded." : "One or more builds failed."
exit(all_ok ? 0 : 1)


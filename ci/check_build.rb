#!/usr/bin/env ruby

require 'fileutils'

#setup log directory
commit_hash = `git rev-parse --short HEAD`.strip
log_dir = "logs/commit-#{commit_hash}"
FileUtils.mkdir_p(log_dir)

#build standalone sim
def build_cpp(label:, build_dir:, compiler: nil, log_file:)
  puts "[#{label}] Building standalone simulator..."

  cmake_cmd = "cmake -S . -B #{build_dir}"
  cmake_cmd += " -DCMAKE_CXX_COMPILER=#{compiler}" if compiler
  cmake_ok = system("#{cmake_cmd} > #{log_file} 2>&1")

  if cmake_ok
    make_ok = system("cmake --build #{build_dir} >> #{log_file} 2>&1")
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
  install_ok = nil
  import_ok = false

  File.open(full_log_path, "w") do |log|
    Dir.chdir(build_dir) do
      IO.popen("pip install .. > /dev/null 2>&1") do |io|
        io.each { |line| log.puts line }
      end
      install_ok = $?.success?
    end
  end

  #find and copy .so into site packages
  built_so = Dir.glob("#{build_dir}/../**/sanafecpp*.so").find { |f| File.file?(f) }
  dest_dir = File.join(ENV["CONDA_PREFIX"], "lib", "python#{RUBY_VERSION[/\d+\.\d+/]}", "site-packages", "sanafe")

  if built_so && File.exist?(built_so)
    FileUtils.mkdir_p(dest_dir)
    FileUtils.cp(built_so, dest_dir)
  end

  #check that import works
  import_output = `python3 -c 'import sanafe' 2>&1`
  import_ok = $?.success?
  File.open(log_file, "a") { |f| f.puts "\n[Import Test]\n#{import_output}" }

  if install_ok && import_ok
    puts "[#{label}] Python build: PASS"
    true
  else
    puts "[#{label}] Python build: FAIL (see #{log_file})"
    false
  end
end

#set clang path
clang_path = "/home/jab23579/useful/clang-20/bin/clang++"

#gcc builds
gcc_sim_ok = build_cpp(
  label: "GCC",
  build_dir: "build_gcc",
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


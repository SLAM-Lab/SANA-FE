#!/usr/bin/env ruby
# Runs the SANA-FE Python unit-test suite (stdlib unittest) and passes only if
# every test completes successfully.
#
# Importing `sanafe` requires the compiled extension to be installed. There is
# no system-wide install: check_build.rb builds and installs the extension into
# per-compiler virtualenvs under the shared CI log directory
# (#{SANAFE_CI_LOG_DIR}/venv_gcc and .../venv_clang). This script reuses those
# interpreters so the API is exercised against the actual built .so -- and, when
# both are present, against both the GCC- and Clang-built extensions.
#
# This check therefore depends on check_build.rb having run earlier in the same
# CI invocation (run.rb sets SANAFE_CI_LOG_DIR so both scripts agree on paths).
require 'fileutils'

log_dir  = ENV["SANAFE_CI_LOG_DIR"] || "logs/commit-latest"
log_file = "#{log_dir}/python.log"
FileUtils.mkdir_p(log_dir)

test_dir  = "tests/python"
test_file = "#{test_dir}/test.py"

# Resolve the interpreters to test against. Prefer the per-compiler venvs that
# check_build.rb creates; these are the only interpreters guaranteed to have
# the freshly built `sanafe` extension installed.
def venv_python(log_dir, label)
  candidate = File.expand_path(File.join(log_dir, "venv_#{label}", "bin", "python"))
  File.exist?(candidate) ? candidate : nil
end

interpreters = []
%w[gcc clang].each do |label|
  py = venv_python(log_dir, label)
  interpreters << [label.upcase, py] if py
end

# Allow an explicit override (e.g. local runs against an interpreter that
# already has sanafe installed, such as a `pip install .` into the active env).
if (override = ENV["SANAFE_PYTHON"])
  interpreters = [["OVERRIDE", override]]
end

def fail_and_exit(log_file, message)
  File.open(log_file, "w") { |log| log.puts message }
  puts message
  exit(1)
end

unless File.exist?(test_file)
  fail_and_exit(log_file,
    "Python Test: FAIL (#{test_file} not found; run from the repository root)")
end

if interpreters.empty?
  fail_and_exit(log_file,
    "Python Test: FAIL (no usable Python interpreter found; expected a build " \
    "venv at #{log_dir}/venv_gcc or #{log_dir}/venv_clang -- has check_build.rb " \
    "run? -- or set SANAFE_PYTHON to an interpreter with sanafe installed)")
end

all_ok = true

File.open(log_file, "w") do |log|
  log.sync = true
  interpreters.each do |label, python|
    # Sanity-check that this interpreter can actually import sanafe before
    # running the suite, so a missing/partial install is reported clearly
    # rather than as a wall of import errors.
    import_ok = system("#{python} -c \"import sanafe\" > /dev/null 2>&1")

    header = "=============================================================\n" \
             "[#{label}] interpreter: #{python}\n" \
             "[#{label}] import sanafe: #{import_ok ? 'OK' : 'FAILED'}\n" \
             "============================================================="
    log.puts header
    puts "[#{label}] Running Python unit tests (#{python})..."

    unless import_ok
      log.puts "[#{label}] Could not 'import sanafe' with this interpreter; " \
               "the extension is not installed in its environment."
      puts "[#{label}] Python Test: FAIL (sanafe not importable; see #{log_file})"
      all_ok = false
      next
    end

    # unittest discovery puts the test directory on sys.path automatically, so
    # the import works regardless of cwd (here, the repo root). A failing or
    # erroring test yields a non-zero exit; skips count as success.
    command = "#{python} -m unittest discover -v -s #{test_dir} -p test.py"
    log.puts "$ #{command}"
    output = `#{command} 2>&1`
    status = $?.exitstatus
    log.puts output

    if status.zero?
      puts "[#{label}] Python Test: PASS"
    else
      puts "[#{label}] Python Test: FAIL (exit #{status}, see #{log_file})"
      all_ok = false
    end
  end

  log.puts "-------------------------------"
  log.puts(all_ok ? "All Python unit tests passed." : "One or more Python test runs failed.")
end

puts all_ok ? "Python Test: PASS" : "Python Test: FAIL (see #{log_file})"
exit(all_ok ? 0 : 1)

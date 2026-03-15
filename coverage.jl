using Coverage, Pkg
# process '*.cov' files
Pkg.test(; coverage = true)
coverage = process_folder()
LCOV.writefile("lcov.info", coverage)
clean_folder("src")

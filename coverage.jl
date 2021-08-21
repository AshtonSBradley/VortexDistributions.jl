using Coverage, Pkg 
# process '*.cov' files
Pkg.test(; coverage=true)
coverage = process_folder() # defaults to src/; alternatively, supply the folder name as argument
LCOV.writefile("lcov.info", coverage)
clean_folder("src")
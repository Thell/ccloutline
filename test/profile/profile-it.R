gpfile <- tempfile(pattern = "ccl-gprof-", tmpdir = ".", fileext = ".out")
Sys.setenv("CPUPROFILE_FREQUENCY" = 250)
Sys.setenv("CPUPROFILE" = gpfile )

library(ccloutline)
cat("\nProfiling...\n")
f <- ccloutline:::.random_source_gen
mat <- f(4096,4096,2)
print(system.time(out <- ccloutline:::ccl_profile(mat)))
cat( "\n",
     paste("google-pprof", "--web --focus=ccl_outlines", "./src/ccloutline.so", gpfile),
     "\n\n" )

library(ccloutline)
library(microbenchmark)

cat("\n")
f <- ccloutline:::.random_source_gen
mat <- f(512,512,2)

tmp <- ccl_labels(mat)[[2]]
cat( "\nTotal connected components: ", length(tmp),
     "\nMembership size summary:\n")
print(summary(tmp))
tmp <- sapply( ccl_edges(mat), function(x) length(x)/2 )
cat( "\nTotal edges: ", sum(tmp), "\nComponent edge summary:\n")
print(summary(tmp))
rm(tmp)

mb <- microbenchmark( borders = ccl_borders(mat),
                      labels = ccl_labels(mat),
                      edges = ccl_edges(mat),
                      cycles = ccl_cycles(mat),
                      paths = ccl_paths(mat),
                      outlines = ccl_outlines(mat))
cat("\n")
print(mb)
mbdiff <- diff( c( 0, summary(mb)$median ) )
cat("\nIsolated units:\n")
names(mbdiff) <- levels(mb$expr)
print( round( rbind( time = mbdiff, percent = mbdiff/sum(mbdiff)), 3) )

# Decisive test (docs/design/ancient-ancestry-rebuild.md §7.10): REAL admixtools2 qpAdm on the same
# EIGENSTRAT (AADR sources + outgroups + huF98AFD target, CHM13, ~20k sites). If WHG resolves to
# ~±2-3%, our pooled-frequency approximation is the bottleneck; if WHG SE is also huge, it's a data
# limit (sparse WHG source).
suppressMessages(library(admixtools))
prefix <- "/Users/jkane/.claude/jobs/f7bfbebc/tmp/qpadm"

left   <- c("WHG","ANF","Steppe")
right  <- c("Han","Papuan","Yoruba","Karitiana","Mbuti","UPEuro","ANE")
target <- "Target"
pops   <- c(left, right, target)

cat("== f2_from_geno (blgsize 0.05 Morgan, maxmiss 0.5) ==\n")
f2 <- f2_from_geno(prefix, pops = pops, maxmiss = 0.5, blgsize = 0.05, verbose = TRUE)

cat("\n== qpadm ==\n")
res <- qpadm(f2, left = left, right = right, target = target)

cat("\n== WEIGHTS ==\n")
print(as.data.frame(res$weights))
cat("\n== MODEL FIT (rankdrop / p-value) ==\n")
print(as.data.frame(res$rankdrop))

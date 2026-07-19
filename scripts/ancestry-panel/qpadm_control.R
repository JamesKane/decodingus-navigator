suppressMessages(library(admixtools))
prefix <- "/Users/jkane/.claude/jobs/f7bfbebc/tmp/aadr446e"
left  <- c("WHG","ANF","Steppe")
right <- c("Han","Papuan","Yoruba","Karitiana","Mbuti","UPEuro","ANE")
f2 <- f2_from_geno(prefix, pops=c(left,right,"English"), maxmiss=0.5, blgsize=0.05, verbose=FALSE)
res <- qpadm(f2, left=left, right=right, target="English")
cat("\n== ENGLISH (1000G British) as target — WEIGHTS ==\n"); print(as.data.frame(res$weights))
cat("\n== MODEL FIT ==\n"); print(as.data.frame(res$rankdrop))

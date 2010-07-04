############################################################
# kmer_hist.r
#
# Produce pdf and txt histograms for further analysis.
############################################################

############################################################
# est.mean
#
# Estimate the coverage mean by finding the max past the
# first valley.
############################################################
est.mean = function(d) {
  # make histogram
  hc = hist(d, breaks=seq(0,round(max(d))+1,1), plot=F)$counts

  # find first valley (or right before it)
  valley = hc[1]
  i = 2
  while(hc[i] < valley) {
    valley = hc[i]
    i = i + 1
  }

  # return max over the rest
  max.hist = hc[i]
  max.cov = i
  while(i <= length(hc)) {
    if(hc[i] > max.hist) {
      max.hist = hc[i]
      max.cov = i
    }
    i = i + 1
  }
  return( max.cov )
}

############################################################
# load data
############################################################
cov = scan("kmers.txt")

cov.est = est.cov(cov)

max.cov = max(cov)

############################################################
# histogram pdf
############################################################
u = est.cov(cov)
pdf("kmer_hist.pdf")
hc = hist(cov, breaks=seq(0,ceiling(max.cov)), plot=F)$counts
hist(cov, breaks=seq(0,ceiling(max.cov)), xlim=c(0, round(1.7*u)), ylim=c(0, round(1.2*hc[u])), xlab="Coverage", main="k-mer counts")
dev.off()

############################################################
# histogram txt
############################################################
scov = cov[cov < 5]
h = hist(scov, breaks=seq(0,5,.1), plot=FALSE)

outf = "kmer_hist.bin10.txt"
cat("", file=outf)
for(i in seq(1,50)) {
  cat(i*.1, "\t", h$counts[i], "\n", sep=" ", file=outf, append=TRUE)
}

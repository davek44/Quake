############################################################
# kmer_model.r
#
# Read in a file full of kmer coverages, and optimize
# a probabilistic model in which kmer coverage is simulated
# as follows:
#
# 1. Choose error or no as binomial(p.e)
# 2. Choose copy number as zeta(p)
# 3a. If error, choose coverage as my discrete pareto(p)
# 3b. If not, choose coverage as normal(copy*u.v, var.v)
#
# Finally, output the ratio of being an error kmer vs being
# a non-error kmer for varying levels of low coverage.
############################################################
library(VGAM)

outf = "cutoff.txt"
max.copy = 30


############################################################
# est.cov
#
# Estimate the coverage mean by finding the max past the
# first valley.
############################################################
est.cov = function(d) {
  hc = d[,2]
  
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
cov = read.table('kmers.hist')[,1:2]
cov.est = est.cov(cov)

# filter extremes from kmers
cov = cov[cov[,1] < 1.25*max.copy*cov.est,]
max.cov = max(cov[,1])


############################################################
# model 
############################################################
model = function(params) {
  print(params)
  
  zp.copy = params[1]
  p.e = params[2]
  shape.e = params[3]
  u.v = params[4]
  var.v = params[5]

  if(zp.copy <= 0 | shape.e <= 0 | var.v <= 0)
    return(list(like=-Inf))
  
  dcopy = dzeta(1:max.copy, p=zp.copy)
  
  kmers.probs = matrix(0, dim(cov)[1], ncol=(max.copy+1))

  # error
  pp = ppareto(1:(max.cov+1), 1, shape.e)
  derr = pp[2:(max.cov+1)] - pp[1:max.cov]
  kmers.probs[,1] = p.e * derr[cov[,1]]

  # no error
  for(copy in 1:max.copy) {
    kmers.probs[,1+copy] = dcopy[copy] * (1-p.e) * dnorm(cov[,1], mean=(copy*u.v), sd=sqrt(copy*var.v))
  }

  logprobs = log(apply(kmers.probs, 1, sum))
  print(sum(logprobs*cov[,2]))
  return(list(like=-sum(logprobs*cov[,2]), probs=kmers.probs))
}

############################################################
# display.params
############################################################
display.params = function(x, print) {
  if(print) {
    cat("zp.copy: ",x[1],"\n", file=outf, append=T)
    cat("p.e: ",x[2],"\n", file=outf, append=T)
    cat("shape.e: ",x[3],"\n", file=outf, append=T)
    cat("u.v: ",x[4],"\n", file=outf, append=T)
    cat("var.v: ",x[5],"\n", file=outf, append=T)
  }
  return(list(zp.copy=x[1], p.e=x[2], shape.e=x[3], u.v=x[4], var.v=x[5]))
}

############################################################
# cutoffs
############################################################
cutoffs = function(p) {
  dcopy = dzeta(1:max.copy, p=p$zp.copy)
  pp = ppareto(1:(max.cov+1), 1, p$shape.e)
  derr = pp[2:(max.cov+1)] - pp[1:max.cov]

  ratios = rep(0,40)  
  for(cov in 1:40) {
    error = p$p.e * derr[cov]

    no.error.list = rep(0,max.copy)
    for(copy in 1:max.copy) {
      no.error.list[copy] = dcopy[copy] * (1-p$p.e) * dnorm(cov, mean=(copy*p$u.v), sd=sqrt(copy*p$var.v))
    }
    
    ratios[cov] = error / sum(no.error.list)
  }
  return(ratios)
}

############################################################
# action
############################################################
init = c(3, .9, 3, cov.est, 5*cov.est)
#ol = c(.001, .001, 0.001, 0, 10)
#ou = c(20, .999, 20, 1000, Inf)
#opt = optim(init, function(x) model(x)$like, lower=ol, upper=ou, method="L-BFGS-B", control=list(trace=1, maxit=1000))
opt = optim(init, function(x) model(x)$like, method="BFGS", control=list(trace=1, maxit=5000))
cat('value:',opt$value,"\n")
p=display.params(opt$par, F)

like.ratio.t = 1000
cut = min((1:40)[cutoffs(p) < like.ratio.t])
cat(cut,"\n", file=outf)

display.params(opt$par, T)
write(cutoffs(p), file=outf, append=T)

max.x = 250
dcopy = dzeta(1:max.copy, p=p$zp.copy)
pp = ppareto(1:(max.cov+1), 1, p$shape.e)
derr = pp[2:(max.cov+1)] - pp[1:max.cov]

kmers.probs = matrix(0, max.x, ncol=(1+max.copy))
# error
kmers.probs[,1] = p$p.e * derr[1:max.x]
# no error
for(copy in 1:max.copy) {
  kmers.probs[,1+copy] = dcopy[copy] * (1-p$p.e) * dnorm(1:max.x, mean=(copy*p$u.v), sd=sqrt(copy*p$var.v))
}

all.dist = apply(kmers.probs, 1, sum)
err.dist = kmers.probs[,1]/p$p.e
norm.dist = apply(kmers.probs[,2:(1+max.copy)], 1, sum)/(1-p$p.e)

ratios = err.dist*p$p.e / norm.dist/(1-p$p.e)

############################################################
# kmer_model_qmer.r
#
# Read in a file full of kmer coverages, and optimize
# a probabilistic model in which kmer coverage is simulated
# as follows:
#
# 1. Choose error or no as binomial(p.e)
# 2. Choose copy number as zeta(zp.copy)
# 3a. If error, choose coverage as gamma(shape.e, scale.e)
# 3b. If not, choose coverage as normal(copy*u.v, copy*var.v)
#
# Finally, output the ratio of being an error kmer vs being
# a non-error kmer for varying levels of low coverage.
#
# Because there is a lot of variation between each integer
# that isn't very gamma-like, I discretize the distribution
# during the optimization.
############################################################
library(VGAM)

# get ratio from input
if(length(commandArgs(trailingOnly=TRUE)) > 0) {
  ratio.goal = as.integer(commandArgs(trailingOnly=TRUE)[1])
} else {
  ratio.goal = 200
}

outf = "cutoff.txt"
max.copy = 30

############################################################
# est.cov
#
# Estimate the coverage mean by finding the max past the
# first valley.
############################################################
est.cov = function(d) {
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

  if(max.hist / valley < 1.5) {
    stop("Cannot identify reliable trusted k-mer coverage peak. Please verify the integrity of your dataset (e.g. by running kmer_hist.r). If there is no clear trusted k-mer distribution, the coverage is too low to correct these reads. If there is, please e-mail dakelley@umiacs.umd.edu and choose the cutoff by hand in the meantime")
  }  
    
  return( max.cov )
}

############################################################
# load data
############################################################
cov = scan("kmers.txt")

cov.est = est.cov(cov)

# filter extremes from kmers
cov = cov[cov < 1.25*max.copy*cov.est]
max.cov = max(cov)
max.cov.discrete = ceiling(max.cov)

# add epsilon to zeroes
cov[cov == 0] = 1e-12

############################################################
# model 
############################################################
model = function(params) {
  #print(params)
  zp.copy = params[1]
  p.e = params[2]
  shape.e = params[3]
  scale.e = params[4]
  u.v = params[5]
  var.v = params[6]

  if(zp.copy <= 0 | p.e <= 0 | p.e >= 1 | shape.e <= 0 | scale.e <= 0 | var.v <= 0)
    return(list(like=-Inf))
  
  dcopy = dzeta(1:max.copy, p=zp.copy)
  
  kmers.probs = matrix(0, nrow=length(cov), ncol=(max.copy+1))

  # error
  dgamma.discrete = pgamma(1:max.cov.discrete, shape=shape.e, scale=scale.e) - pgamma(0:(max.cov.discrete-1), shape=shape.e, scale=scale.e)
  kmers.probs[,1] = p.e * dgamma.discrete[ceiling(cov)]
  #kmers.probs[,1] = p.e * dgamma(cov, shape=shape.e, scale=scale.e)

  # no error
  for(copy in 1:max.copy) {
    kmers.probs[,1+copy] = dcopy[copy] * (1-p.e) * dnorm(cov, mean=(copy*u.v), sd=sqrt(copy*var.v))
  }

  logprobs = log(apply(kmers.probs, 1, sum))
  #print(sum(logprobs))
  return(list(like=-sum(logprobs), probs=kmers.probs))
}


############################################################
# display.params
############################################################
display.params = function(x, print) {
  if(print) {
    cat("zp.copy: ",x[1],"\n", file=outf, append=T)
    cat("p.e: ",x[2],"\n", file=outf, append=T)
    cat("shape.e: ",x[3],"\n", file=outf, append=T)
    cat("scale.e: ",x[4],"\n", file=outf, append=T)
    cat("u.v: ",x[5],"\n", file=outf, append=T)
    cat("var.v: ",x[6],"\n", file=outf, append=T)
  }
  return(list(zp.copy=x[1], p.e=x[2], shape.e=x[3], scale.e=x[4], u.v=x[5], var.v=x[6]))
}

############################################################
# cutoff
#
# Use a grid search to find the coverage number that
# makes the ratio of error probability to nonerror
# probability closest to ratio.goal
############################################################
ratios = function(cov, p) {  
  dcopy = dzeta(1:max.copy, p=p$zp.copy)
  error = p$p.e * dgamma(cov, shape=p$shape.e, scale=p$scale.e)
  no.error = 0
  for(copy in 1:max.copy) {
    no.error = no.error + dcopy[copy] * (1-p$p.e) * dnorm(cov, mean=(copy*p$u.v), sd=sqrt(copy*p$var.v))
  }
  return(error/no.error)
}

error.max = function(p) {
  cov.max = 0
  dens.max = 0
  for(c in seq(0,30, .01)) {
    dens = dgamma(c, shape=p$shape.e, scale=p$scale.e)
    if(dens > dens.max) {
      dens.max = dens
      cov.max = c
    }
  }
  return(cov.max)
}

cutoff = function(p) {
  cut.start = error.max(p)
  
  #ratio.goal = exp(3 + .05*p$u.v)
  #ratio.goal = 1000
  cov.best = 0
  ratio.best = abs(ratios(0,p) - ratio.goal)
  for(c in seq(cut.start, 25, .01)) {
    r = abs(ratios(c,p) - ratio.goal)
    if(r < ratio.best) {
      ratio.best = r
      cov.best = c
    }
    #cat(c," ",ratios(c,p),"\n")
  }
  return(cov.best)
}

############################################################
# action
############################################################
init = c(2, .9, .8, .5, cov.est, 3*cov.est)
#ol = c(.001, .001, 0.001, 0, 10)
#ou = c(20, .999, 20, 1000, Inf)
#opt = optim(init, function(x) model(x)$like, lower=ol, upper=ou, method="L-BFGS-B", control=list(trace=1, maxit=1000))
opt = optim(init, function(x) model(x)$like, method="BFGS", control=list(trace=1, maxit=1000))
cat('value:',opt$value,"\n")
p=display.params(opt$par, F)

cut = cutoff(p)
cat(cut,"\n", file=outf)
display.params(opt$par, T)
write(cut, file=outf, append=T)

max.x = 250
dcopy = dzeta(1:max.copy, p=p$zp.copy)
  
kmers.probs = matrix(0, nrow=max.x, ncol=(max.copy+1))
# error
#kmers.probs[,1] = p$p.e * (ppareto(2:251, loc.e, p$shape.e) - ppareto(1:250, loc.e, p$shape.e))
kmers.probs[,1] = p$p.e * (pgamma(1:250, shape=p$shape.e, scale=p$scale.e) - pgamma(0:249, shape=p$shape.e, scale=p$scale.e))

# no error
for(copy in 1:max.copy) {
  kmers.probs[,1+copy] = dcopy[copy] * (1-p$p.e) * (pnorm(1:250, mean=(copy*p$u.v), sd=sqrt(copy*p$var.v)) - pnorm(0:249, mean=(copy*p$u.v), sd=sqrt(copy*p$var.v)))
}

all.dist = apply(kmers.probs, 1, sum)
err.dist = kmers.probs[,1]/p$p.e
norm.dist = apply(kmers.probs[,2:(1+max.copy)], 1, sum)/(1-p$p.e)

en.ratios = err.dist*p$p.e/norm.dist/(1-p$p.e)

####################
## Bart KM Jacobs ##
####################

rundPCR <- function(lambda,mudropc=20000,sddropc=0,mudropr=1,sddropr=0,Pvar=T,piperr=0,seed=runif(1),dropsd=0,falpos=0,falneg=0,reps=8,sims=1000,fullres=F){
# lambda: vector or number of the true concentration (as expected target copies/partition)
# mudropc: average number of partitions generated
#		Default 20000 based on the Bio-Rad ddPCR QX100 and Life QuantStudio expected values
# sddropc: standard deviation of the number of partitions generated
#		Default 0 for constant number of partitions
# mudropr: average proportion (between 0 and 1) of retained partitions
#		Default 1 for no loss (on average)
# sddropr: standard deviation of the proportion of retained partitions
#		Default 0 for a constant loss
### It is assumed that both generating and retaining of partitions is ~N
# Pvar: TRUE: number of copies in constant volume follows P(c) distribution
#	      FALSE: number of copies in constant volume is constant
#   Default TRUE for the realistic Poisson model.
# piperr: coefficient of variation of the actually pipetted volume.
#		Default 0 for constant volume equal to the expected volume
# seed: random seed to obtain the same results if rerun
# dropsd: relative variability of the partition volume
#         parameter sigma of a lognormal distribution with mu=0
#		Default 0 for constant partition size
# falpos: probability that a partition containing no copy gives a positive result
#		Default 0 for no false positives
# falneg: probability that a partition containing at least one copy gives a negative result
#		Default 0 for no false negatives
# repw: number of technical replicates in each run.
#   Default 8 based on classic 12x8 plates
# sims: number of simulations for each given lambda
#   Default 1000
# fullres: output stored, statistics for every similation (TRUE) or only summary statistics (FALSE)
#   Default FALSE for faster speed
  
set.seed(seed)
if(max(!is.finite(lambda))) stop("Concentrations should all be numeric")
if(min(lambda)<0) stop("lambda cannot be negative")

# Same procedure for all replicates
# repfunc is called internally in samfunc
repfunc <- function(repdat){
  dropmem <- sample(repdat[1],repdat[2],replace=T,prob=rlnorm(repdat[1],0,dropsd))
    # partition membership, probability proportional to size, size following a lognormal distribution
  dropn <- ifelse(mudropr>=1,repdat[1],round(repdat[1]*plogis(rnorm(1,log(mudropr/(1-mudropr)),log((mudropr+sddropr)/(mudropr-sddropr)*(1-mudropr+sddropr)/(1-mudropr-sddropr))/2))))
    # number of partitions retained
  dropmem <- dropmem[dropmem<=dropn]
    # only retain copies of which the partition is retained (lower rank)
  dropno <- dropn - length(as.vector(table(dropmem)))
    # number of partitions without copy (total - number with copies)
  dropneg <- rbinom(1,dropno,1-falpos)+rbinom(1,dropn-dropno,falneg)
    # number of partitions with a negative signal (true neg + false neg)
  lambdahat <- -log(dropneg/dropn)
    # estimate per replicate
  sd.asym <- sqrt(1/dropneg-1/dropn)
    # asymptotic standard deviation
  BIlow <- max(lambdahat-1.96*sd.asym,0)
    # lower bound confidence interval
  BIup <- lambdahat+1.96*sd.asym
    # upper bound confidence interval
  BIw <- BIup-BIlow
    # width of confidence interval
  inBI <- (repdat[3] > BIlow)*(repdat[3] < BIup)
    # 0/1 indicator 1 true value in CI, 0 not.
  returnv <- c(dropn,dropneg,lambdahat,sd.asym,BIw,inBI)
  return(returnv)
    # returnv is a vector with the following elements:
    # 1. analyzed partitions 2. negative partitions 3. concentration estimate
    # 4. asymptotic standard deviation 5. CI width 6. in CI indicator
}

        
# Same procedure for all simulations
# samfunc is called internally in lamfunc
samfunc <- function(lambdan){
      dropstart <- round(rnorm(reps,mudropc,sddropc))
        # number of partitions
      copyvar <- lambdan*rnorm(reps,1,piperr)*dropstart
	copyvar[copyvar<0] <- 0
        # expected number of copies after pipette variation
	  	copyn <- ifelse(rep(Pvar,reps),rpois(reps,copyvar),round(copyvar))
        # number of copies
      lamdummy <- rep(lambdan,reps)
      repdat <- as.list(data.frame(rbind(dropstart,copyn,lamdummy)))
        # number of partitions and copies in a list with reps elements, all pairs
        # partitions is the first element, copies the second
  repres <- sapply(repdat,repfunc)
  droptot <- sum(repres[1,])
    # total number of partitions over all 8 replicates
  dropnegtot <- sum(repres[2,])
    # total number of negative partitions over all 8 replicates
  M1.lambdahat <- -log(dropnegtot/droptot)
    # pooled estimate
  M1.sd.asym <- sqrt(1/dropnegtot-1/droptot)
    # pooled asymptotic standard deviation
  M1.BIlow <- max(M1.lambdahat-1.96*M1.sd.asym,0)
    # lower bound confidence interval
  M1.BIup <- M1.lambdahat+1.96*M1.sd.asym
    # upper bound confidence interval
  M1.BIw <- M1.BIup-M1.BIlow
    # width of confidence interval
  M1.inBI <- (lambdan > M1.BIlow)*(lambdan < M1.BIup)
    # 0/1 indicator 1 true value in CI, 0 not.
  M2.lambdahat <- mean(repres[3,])
    # replicate based estimate
  M2.sd <- sd(repres[3,])
    # emperical standard deviation
  M2.BIlow <- max(M1.lambdahat-qt(0.975,reps-1)*M2.sd/sqrt(reps),0)
    # lower bound confidence interval
  M2.BIup <- M2.lambdahat+qt(0.975,reps-1)*M2.sd/sqrt(reps)
    # upper bound confidence interval
  M2.BIw <- M2.BIup-M2.BIlow
    # width of confidence interval
  M2.inBI <- (lambdan > M2.BIlow)*(lambdan < M2.BIup)
    # 0/1 indicator 1 true value in CI, 0 not.
  if(fullres==T){
    returnv <- c(M1.lambdahat,M1.sd.asym,M1.BIw,M1.inBI,M2.lambdahat,M2.sd,M2.BIw,M2.inBI,t(repres[3:6,]),t(repres[1:2,]),droptot,dropnegtot)
      # Vector with 6*reps+10 elements
      # 1. - 8.: see fullres==F
      # 9. - reps+8: lambdahat for replicates
      # reps+9 - 2*reps+8: sd.asym for replicates
      # 2*reps+9 - 3*reps+8: CIwidth replicates
      # 3*reps+9 - 4*reps+8: 1/0 indicators in CI replicates
      # 4*reps+9 - 5*reps+8: number of partitions
      # 5*reps+9 - 6*reps+8: number of negative partitions
      # 6*reps+9: pooled total number of partitions
      # 6*reps+10: pooled total number of negative partitions
  }
  if(fullres==F){
    returnv <- c(M1.lambdahat,M1.sd.asym,M1.BIw,M1.inBI,M2.lambdahat,M2.sd,M2.BIw,M2.inBI)
      # Vector with 8 elements:
      # 1. lambdahat pooled, 2. sd pooled, 3. CI width pooled, 4. in CI pooled
      # 5. lambdahat replicate based, 6. sd replicate based, 7. CI width replicate based, 8. in CI replicate based
  }
  return(returnv)
}
  

# Function for the lambda's
lamfunc <- function(lambdai){
  simvec <- rep(lambdai,sims)
  simres <- sapply(simvec,samfunc)
  if(fullres==T){
    returnv <- t(simres)
    # Vector with (6*reps+10)*sims elements
    # 1. - sims: lambdahat pooled
    # sims+1 - 2*sims: sd pooled
    # 2*sims+1 - 3*sims: CI width pooled
    # 3*sims+1 - 4*sims: in CI pooled
    # 4*sims+1 - 5*sims: lambdahat replicate based
    # 5*sims+1 - 6*sims: sd replicate based
    # 6*sims+1 - 7*sims: CI width replicate based
    # 7*sims+1 - 8*sims: in CI replicate based
    # 8*sims+1 - (reps+8)*sims: lambdahat for replicates
    # (reps+8)*sims+1 - (2*reps+8)*sims: sd.asym for replicates
    # (2*reps+8)*sims+1 - (3*reps+8)*sims: CIwidth replicates
    # (3*reps+8)*sims+1 - (4*reps+8)*sims: 1/0 indicators in CI replicates
    # (4*reps+8)*sims+1 - (5*reps+8)*sims: number of partitions
    # (5*reps+8)*sims+1 - (6*reps+8)*sims: number of negative partitions
    # (6*reps+8)*sims+1 - (6*reps+9)*sims: pooled total number of partitions
    # (6*reps+9)*sims+1 - (6*reps+10)*sims: pooled total number of negative partitions
  }
  if(fullres==F){
    M1.lambdahat <- mean(simres[1,])
      # 1. average pooled estimate
    M1.simsd <- sd(simres[1,])
      # 2. standard deviation among pooled estimates
    M1.essd <- mean(simres[2,])
      # 3. estimated standard deviation of a pooled estimate
    M1.cover <- mean(simres[4,])
      # 4. coverage pooled estimates
    M2.lambdahat <- mean(simres[5,])
      # 5. average replicate based estimate
    M2.simsd <- sd(simres[5,])
      # 6. standard deviation among replicate based estimates
    M2.essd <- mean(simres[6,])/sqrt(8)
      # 7. estimated standard deviation of a replicate based estimate
    M2.cover <- mean(simres[8,])
      # 8. coverage replicate based estimates    
    returnv <- c(M1.lambdahat,M1.simsd,M1.essd,M1.cover,M2.lambdahat,M2.simsd,M2.essd,M2.cover)
      # Vector with 8 elements described above.
  }
  return(returnv)
}


out <- sapply(lambda,lamfunc)


if(fullres==T){
  lambdahat <- out[1:sims,]
    # estimates for all simulations
  sd.asym <- out[(sims+1):(2*sims),]
    # pooled asymptotic standard deviations for all simulations
  CI.width <- out[(2*sims+1):(3*sims),]
    # confidence interval width for all simulations
  in.CI <- out[(3*sims+1):(4*sims),]
    # 0/1 indicator whether true parameter in the interval for all simulations
  parts <- out[((6*reps+8)*sims+1):((6*reps+9)*sims),]
    # total number of analyzed partitions for every replicates in all simulations
  nparts <- out[((6*reps+9)*sims+1):((6*reps+10)*sims),]
    # number of analyzed partitions with a negative signal for every replicates in all simulations
  estimate <- colMeans(lambdahat)
    # average of the pooled estimates over all simulations (=lambda?)
  bias <- estimate - lambda
    # bias of the pooled estimate
  sd <- (colSums(lambdahat^2) - colSums(lambdahat)^2/sims) /(sims-1)
    # standard deviation among the pooled estimates (NOT the asymptotic standard deviation!)
  RMSE <- sqrt(bias^2 + sd^2)
    # root mean squared error
  bias.stdz <- bias/lambda
    # standardized bias of the pooled estimate
  sd.stdz <- sd/lambda
    # standardized standard deviation of the pooled estimate
  RMSE.stdz <- RMSE/lambda
    # standardized RMSE of the pooled estimate
  coverage <- colMeans(in.CI)
    # coverage of the pooled estimate
  pooled <- list(estimate=estimate,bias=bias,sd=sd,RMSE=RMSE,coverage=coverage,bias.stdz=bias.stdz,sd.stdz=sd.stdz,RMSE.stdz=RMSE.stdz,lambdahat=lambdahat,sd.asym=sd.asym,CI.width=CI.width,in.CI=in.CI,parts=parts,nparts=nparts)
    # list with all results of the pooled estimation procedure

  lambdahat <- out[(4*sims+1):(5*sims),]
    # estimates for all simulations
  sd.emp <- out[(5*sims+1):(6*sims),]
    #  emperical standard deviations for all simulations
  CI.width <- out[(6*sims+1):(7*sims),]
    # confidence interval width for all simulations
  in.CI <- out[(7*sims+1):(8*sims),]
    # 0/1 indicator whether true parameter in the interval for all simulations
  estimate <- colMeans(lambdahat)
    # average of the replicate based estimates over all simulations (=lambda?)
  bias <- estimate - lambda
    # bias of the replicate based estimate
  sd <- (colSums(lambdahat^2) - colSums(lambdahat)^2/sims) /(sims-1)
    # standard deviation among the replicate based estimates
  RMSE <- sqrt(bias^2 + sd^2)
    # root mean squared error
  bias.stdz <- bias/lambda
    # standardized bias of the replicate based estimate
  sd.stdz <- sd/lambda
    # standardized standard deviation of the replicate based estimate
  RMSE.stdz <- RMSE/lambda
    # standardized RMSE of the replicate based estimate
  coverage <- colMeans(in.CI)
    # coverage of the replicate based estimate
  repbased <- list(estimate=estimate,bias=bias,sd=sd,RMSE=RMSE,coverage=coverage,bias.stdz=bias.stdz,sd.stdz=sd.stdz,RMSE.stdz=RMSE.stdz,lambdahat=lambdahat,sd.emp=sd.emp,CI.width=CI.width,in.CI=in.CI)
    # list with all results of the replicate based estimation procedure
  
  lambdahat <- out[(8*sims+1):((reps+8)*sims),]
    # estimates for all replicates in all simulations
  sd.asym <- out[((reps+8)*sims+1):((2*reps+8)*sims),]
    # asymptotic standard deviations for all replicates in all simulations
  CI.width <- out[((2*reps+8)*sims+1):((3*reps+8)*sims),]
    # confidence interval width for all replicates in all simulations
  in.CI <- out[((3*reps+8)*sims+1):((4*reps+8)*sims),]
    # 0/1 indicator whether true parameter in the interval for all replicates in all simulations
  parts <- out[((4*reps+8)*sims+1):((5*reps+8)*sims),]
    # total number of analyzed partitions for every replicates in all simulations
  nparts <- out[((5*reps+8)*sims+1):((6*reps+8)*sims),]
    # number of analyzed partitions with a negative signal for every replicates in all simulations
  estimate <- colMeans(lambdahat)
    # average of the estimates over all replicates in all simulations (=lambda?)
  bias <- estimate - lambda
    # bias of the estimate
  sd <- (colSums(lambdahat^2) - colSums(lambdahat)^2/(reps*sims)) /(reps*sims-1)
    # standard deviation among the estimates (NOT the asymptotic standard deviation!)
  RMSE <- sqrt(bias^2 + sd^2)
    # root mean squared error
  bias.stdz <- bias/lambda
    # standardized bias of the estimate
  sd.stdz <- sd/lambda
    # standardized standard deviation of the estimate
  RMSE.stdz <- RMSE/lambda
    # standardized RMSE of the estimate
  coverage <- colMeans(in.CI)
    # coverage of the estimate
  singlesample <- list(estimate=estimate,bias=bias,sd=sd,RMSE=RMSE,coverage=coverage,bias.stdz=bias.stdz,sd.stdz=sd.stdz,RMSE.stdz=RMSE.stdz,lambdahat=lambdahat,sd.asym=sd.asym,CI.width=CI.width,in.CI=in.CI,parts=parts,nparts=nparts)
    # list with all results of the estimation procedure for all replicates in all simulations
  
  output <- list(repbased=repbased,pooled=pooled,singlesample=singlesample)
}

if(fullres==F){
  estimate <- out[1,]
    # average of the pooled estimates over all simulations (estimate for lambda)
  bias <- estimate - lambda
    # bias of the pooled estimate
  sd <- out[2,]
    # standard deviation among the pooled estimates (NOT the asymptotic standard deviation!)
  RMSE <- sqrt(bias^2 + sd^2)
    # root mean squared error
  essd <- out[3,]
    # average estimated standard deviation of the pooled estimate
  bias.stdz <- bias/lambda
    # standardized bias of the pooled estimate
  sd.stdz <- sd/lambda
    # standardized standard deviation of the pooled estimate
  RMSE.stdz <- RMSE/lambda
    # standardized RMSE of the pooled estimate
  essd.stdz <- essd/lambda
    # standardized average estimated standard deviation of the pooled estimate
  coverage <- out[4,]
    # coverage of the pooled estimate
  pooled <- list(estimate=estimate,bias=bias,sd=sd,RMSE=RMSE,essd=essd,coverage=coverage,bias.stdz=bias.stdz,sd.stdz=sd.stdz,RMSE.stdz=RMSE.stdz,essd.stdz=essd.stdz)
    # list with summarized results of the pooled estimation procedure
  
  estimate <- out[5,]
    # average of the replicate based estimates over all simulations (estimate for lambda)
  bias <- estimate - lambda
    # bias of the replicate based estimate
  sd <- out[6,]
    # standard deviation among the replicate based estimates
  RMSE <- sqrt(bias^2 + sd^2)
    # root mean squared error
  essd <- out[7,]
    # average estimated standard deviation of the replicate based estimate
  bias.stdz <- bias/lambda
    # standardized bias of the replicate based estimate
  sd.stdz <- sd/lambda
    # standardized standard deviation of the replicate based estimate
  RMSE.stdz <- RMSE/lambda
    # standardized RMSE of the replicate based estimate
  essd.stdz <- essd/lambda
    # standardized average estimated standard deviation of the replicate based estimate
  coverage <- out[8,]
    # coverage of the replicate based estimate
  repbased <- list(estimate=estimate,bias=bias,sd=sd,RMSE=RMSE,essd=essd,coverage=coverage,bias.stdz=bias.stdz,sd.stdz=sd.stdz,RMSE.stdz=RMSE.stdz,essd.stdz=essd.stdz)
    # list with all results of the replicate based estimation procedure
  
  output <- list(repbased=repbased,pooled=pooled)
}

  par(mfrow=c(2,2))
  par(mar=c(4,4,1,1))
  try(plot(lambda,repbased$sd.stdz,log="xy",ylim=c(min(repbased$sd.stdz),max(repbased$sd.stdz)),xlab=expression(lambda),ylab="Standardized standard deviation"))
  try(plot(lambda,repbased$bias.stdz+1,log="xy",ylim=c(min(repbased$bias.stdz+1),max(repbased$bias.stdz+1)),xlab=expression(lambda),ylab="Standardized bias"))
  try(abline(h=1))
  try(plot(lambda,repbased$RMSE.stdz,log="xy",ylim=c(min(repbased$RMSE.stdz),max(repbased$RMSE.stdz)),xlab=expression(lambda),ylab="Standardized RMSE"))
  try(plot(lambda,repbased$coverage,log="x",ylim=c(0,1),xlab=expression(lambda),ylab="95% CI Coverage"))
  try(abline(h=0.95))
  par(mfrow=c(1,1))

return(output)
}


npoptimal <- function(r,theta,n,B,rseed){
    
    set.seed(rseed) 
    select=rep(0,length(r))
    j=1
    while(j<=B){
      latent=runif(n,0,1)
      tox=matrix(nrow=n,ncol=length(r))
      i=1
      while(i<=n){
        tox[i,]<-as.numeric(latent[i]<=r)
        i=i+1
      }
      loss=abs(colMeans(tox)-theta)
      sugglev=which(loss==min(loss))
      mtd=ifelse(length(sugglev)==1,sugglev,sample(sugglev,1))
      select[mtd]=select[mtd]+1
      j=j+1
    }
    cat("True DLT probability:     ", format(round(r,3), nsmall = 2), sep="\t",  "\n");
    cat("MTD selection percentage: ", formatC((select/B)*100, digits=1, format="f"), sep="\t",  "\n");
    cat("Accuracy Index:           ", round(1-length(r)*(sum(abs(r-theta)*(select/B))/sum(abs(r-theta))),4), sep="\t",  "\n");
  }
  

##True toxicity scenarios
r1=c(0.05,0.15,0.30,0.40,0.50,0.60)
r2=c(0.08,0.12,0.20,0.30,0.42,0.53)
r3=c(0.30,0.38,0.45,0.55,0.70,0.80)
r4=c(0.02,0.05,0.10,0.15,0.23,0.30)
r5=c(0.18,0.28,0.36,0.44,0.52,0.65)
r6=c(0.01,0.03,0.05,0.12,0.30,0.46)


####################################
#
#	Input:
#
#	r = true toxicity probabilities 
#	n = sample size
#	theta = target toxicity rate
#	B = number of simulated trials
#	rseed = seed of the random number generator
# 
####################################

theta=0.30
n=25
B=10000
r=r4
rseed=352780

npoptimal(r,theta,n,B,rseed)

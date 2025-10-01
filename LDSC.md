# Load in the packages
```{r}
#install.packages("devtools")
library(devtools)
#install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
```

### Run multivariable LDSC
```{r}
#vector of munged summary statisitcs
traits<-c("ADHD.sumstats.gz","VD.sumstats.gz", "BRH.sumstats.gz", "Overall.sumstats.gz", "MA.sumstats.gz", "MO.sumstats.gz")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(.5,.5,.5,.5,.5,.5)

#vector of population prevalences
population.prev<-c(.05,.15,.15,.15,.1,.2)

#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<-c("ADHD","VD", "BRH", "Overall", "MA", "MO")

#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)

#optional command to save the output as a .RData file for later use
save(LDSCoutput,file="LDSCoutput.RData")
```

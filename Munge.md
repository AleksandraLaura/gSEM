# Load in the packages
```{r}
#install.packages("devtools")
library(devtools)
#install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
```

### Munge data
```{r}
#create vector of the summary statistics files
files<-c("ADHD_withrsID.txt", "Visual_disturbance_hg37.txt", "Bad_recurrent_headaches_hg37.txt", "Migrain_all_hg37.txt", "Migraine_with_aura_hg37.txt", "Migraine_without_aura_hg37.txt")

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
hm3<-"eur_w_ld_chr/w_hm3.snplist"

#name the traits 
trait.names<-c("ADHD","VD", "BRH", "Overall", "MA", "MO")

#list the sample sizes.
N=c(NA, 89653.16, 146060.64, 299106.14, 65597.11, 46462.87)

#definte the imputation quality filter
info.filter=0.9

#define the MAF filter
maf.filter=0.01

#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)
```

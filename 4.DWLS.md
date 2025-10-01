# Load in the packages
```{r}
#install.packages("devtools")
library(devtools)
#install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
```

### Run using DWLS
```{r}
#To run using DWLS estimation
CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")

#print CommonFactor_DWLs output
CommonFactor_DWLS
```

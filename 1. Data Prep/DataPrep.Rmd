# Load in the packages
```{r}
#install.packages("devtools")
library(devtools)
#install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
```

### Lift over migraine data from hg38 to hg19
```{r}
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("rtracklayer")

library(rtracklayer)
library(dplyr)
```

#### Did this for each of the five files

##### VD
```{r}
# read your table
df <- read.table("Visual_disturbance.cl.txt", header = TRUE, stringsAsFactors = FALSE)

# convert to GRanges (0-based start)
gr <- GRanges(
  seqnames = df$Chr,
  ranges = IRanges(start = df$PosB38, end = df$PosB38),
  names = df$ID
)

# load chain file
chain <- import.chain("hg38ToHg19.over.chain")

# liftover
gr37 <- liftOver(gr, chain)

# add hg19 positions to your table
df$Pos <- start(gr37)
df$Chr <- as.character(seqnames(gr37))
df = subset(df, select = -c(PosB38, ID) )
df$Chr <- sub("^chr", "", df$Chr)

# write out updated table
write.table(df, "Visual_disturbance_hg37.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

##### BRH
```{r}
# read your table
df <- read.table("Bad_recurrent_headaches.cl.txt", header = TRUE, stringsAsFactors = FALSE)

# convert to GRanges (0-based start)
gr <- GRanges(
  seqnames = df$Chr,
  ranges = IRanges(start = df$PosB38, end = df$PosB38),
  names = df$ID
)

# load chain file
chain <- import.chain("hg38ToHg19.over.chain")

# liftover
gr37 <- liftOver(gr, chain)

# add hg19 positions to your table
df$Pos <- start(gr37)
df$Chr <- as.character(seqnames(gr37))
df = subset(df, select = -c(PosB38, ID) )
df$Chr <- sub("^chr", "", df$Chr)

# write out updated table
write.table(df, "Bad_recurrent_headaches_hg37.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

##### Overall
```{r}
# read your table
df <- read.table("Migrain_all.cl.txt", header = TRUE, stringsAsFactors = FALSE)

# convert to GRanges (0-based start)
gr <- GRanges(
  seqnames = df$Chr,
  ranges = IRanges(start = df$PosB38, end = df$PosB38),
  names = df$ID
)

# load chain file
chain <- import.chain("hg38ToHg19.over.chain")

# liftover
gr37 <- liftOver(gr, chain)

# add hg19 positions to your table
df$Pos <- start(gr37)
df$Chr <- as.character(seqnames(gr37))
df = subset(df, select = -c(PosB38, ID) )
df$Chr <- sub("^chr", "", df$Chr)

# write out updated table
write.table(df, "Migrain_all_hg37.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

##### MA
```{r}
# read your table
df <- read.table("Migraine_with_aura.cl.txt", header = TRUE, stringsAsFactors = FALSE)

# convert to GRanges (0-based start)
gr <- GRanges(
  seqnames = df$Chr,
  ranges = IRanges(start = df$PosB38, end = df$PosB38),
  names = df$ID
)

# load chain file
chain <- import.chain("hg38ToHg19.over.chain")

# liftover
gr37 <- liftOver(gr, chain)

# add hg19 positions to your table
df$Pos <- start(gr37)
df$Chr <- as.character(seqnames(gr37))
df = subset(df, select = -c(PosB38, ID) )
df$Chr <- sub("^chr", "", df$Chr)

# write out updated table
write.table(df, "Migraine_with_aura_hg37.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

##### MO
```{r}
# read your table
df <- read.table("Migraine_without_aura.cl.txt", header = TRUE, stringsAsFactors = FALSE)

# convert to GRanges (0-based start)
gr <- GRanges(
  seqnames = df$Chr,
  ranges = IRanges(start = df$PosB38, end = df$PosB38),
  names = df$ID
)

# load chain file
chain <- import.chain("hg38ToHg19.over.chain")

# liftover
gr37 <- liftOver(gr, chain)

# add hg19 positions to your table
df$Pos <- start(gr37)
df$Chr <- as.character(seqnames(gr37))
df = subset(df, select = -c(PosB38, ID) )
df$Chr <- sub("^chr", "", df$Chr)

# write out updated table
write.table(df, "Migraine_without_aura_hg37.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

### Clean AHDH data
```{r}
ADHD<-read.table("daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta/daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta", header = TRUE)
```

```{r}
#load in stringr package for using strsplit function below
require(stringr)


#reformat SNP column to only include rsID
ADHD$SNP <- sub("^chr", "", sapply(strsplit(ADHD$SNP, ":"), `[`, 1))

#output datafile
write.table(ADHD, file = "ADHD_withrsID.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)
```


### Get Migraine Cohorts
```{r}
Migraine_cohorts <- read.csv("migraine_cohorts.csv", header = TRUE)

Migraine_cohorts$v <- Migraine_cohorts$cases / (Migraine_cohorts$cases + Migraine_cohorts$controls)

Migraine_cohorts$EffN <- 4*Migraine_cohorts$v*(1-Migraine_cohorts$v)*(Migraine_cohorts$cases + Migraine_cohorts$controls)

sum(Migraine_cohorts$EffN)
```




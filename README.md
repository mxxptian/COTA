# COTA

`COTA` is an R package to identify candidate core disease genes using trans-regulatory effects.

# Installation


You can install the development version of
`COTA` from Github via the `devtools` package. I suppose using
the `remotes` package would work as well.

Before installation of COTA, you are also requested the below packages:
``` r
install.packages(c('qvalue', 'data.table', 'stringr', 'tidyr', 'AnnotationDbi', 'org.Hs.eg.db', 'ggplot2', 'igraph', 'VennDiagram', 'biomaRt', 'plyr', 'dplyr'), dependencies=TRUE)

```

``` r
devtools::install_github("mxxptian/COTA")
```

# Example

You can access the example code and data through: https://www.dropbox.com/scl/fo/b1r27fxqlq84ywory8qmo/ABRuFJpz1FIAy9vluMhLAMU?rlkey=n3pzgkdd0fqz9waamidlyfcfy&dl=0
```r
library(qvalue)
library(data.table)
library(stringr)
library(plyr)
library(tidyr)
library(tidyverse)
library(readr)
library(gprofiler2)
library(VennDiagram)
library(igraph)
library(COTA)
library(biomaRt)
library(rsnps)



################################# Gene->Gene->Trait CASE ######################################
# Here we use UKBB summary statistics and GBAT result as an example for gene-gene case

# load UKBB summary statistics
dat = read.table(gzfile('GCST90082964_buildGRCh38.tsv.gz'), header = TRUE)


dat_M3.1 = dat[dat$effect_allele=='M3.1',]



head(dat_M3.1)

out.pheno <- strsplit(as.character(dat_M3.1$Name),'\\.') 

info.out.pheno = do.call(rbind, out.pheno)

out.new <- strsplit(as.character(info.out.pheno[,1]),'\\(')

info.out.new = do.call(rbind, out.new)

dat_M3.1$Name = info.out.new[,1]

dat_M3.1_new = dat_M3.1[!dat_M3.1$Name%in%names(which(table(dat_M3.1$Name)!=1)),]


head(dat_M3.1_new)

cau = data.frame(pval = dat_M3.1_new$p_value)
rownames(cau) = dat_M3.1_new$Name


p.wgs <- cau[,1]
p.wgs <- na.omit(p.wgs)
names(p.wgs) = dat_M3.1_new$Name




# load GBAT data

p.trans <- read.delim('/GBAT/trans_pval/cor_chr1_all.txt', sep = '\t')

for (i in 2:22) {
  temp = read.delim(paste0('/GBAT/trans_pval/cor_chr',i,'_all.txt'), sep = '\t')
  
  p.trans = cbind(p.trans,temp)
}
rm(temp)



threshold = 0.05/nrow(cau)
eta.wgs = threshold

# load reference table including gene position and gene type. Note!!! this position is GRCH 37


ref.table = read.delim('gene_position.txt', sep = '\t')

## only keep gene type in (lincRNA, protein_coding)

ref.table.keep <- ref.table[ref.table$type %in% c('lincRNA', 'protein_coding'),]

## remove genes on chr M, X, and Y
ref.table.keep <- ref.table.keep[!(ref.table.keep$Chromosome %in% c('chrM', 'chrX', 'chrY')),]

ref.table.keep <-ref.table.keep[!duplicated(ref.table.keep$gene_name),]


## standardize gene name
out.ref <- strsplit(as.character(ref.table.keep$gene_name),'\\.')

out.ref.new = do.call(rbind, out.ref)

ref.table.keep$gene_name = out.ref.new[,1]

# Here, we consider all gene1 included in UKBB. You may specify a vector of candidate genes 1

candidate = dput(colnames(p.trans))

# STEP 1: we applied COTA with UKBB and GBAT data with target fdr level at 0.1
result = med_gene(p.trans, p.wgs, ref.table, candidate, target.fdr=0.1, dist=5e6, gene1.type = 'Gene')




# STEP 2: we obatined the pair identified by COTA
result.pair = calc_pair.gene(result$mat.sig, result$mat.p,p.wgs, result$gene1, ref.table.keep,
                                     eta.wgs=threshold)


# STEP 3: we plotted the pair identified by COTA
setwd('/Users/px/Desktop/test/example/')

## load conservation score table
conv.table = read_excel('conservation_scores.xlsx')

## set saving directory of the figures
pic_dir = '/Users/px/Desktop/test/example/'



gen_fig(result.pair$gene.pair, result.pair$sig_gene2, result.pair$non_sig.gene2, p.wgs, eta.wgs=1e-5, conv.table, pic_dir)




################################# SNP->Gene->Trait CASE ######################################



# Here we use UKBB summary statistics and GBAT result as an example for snp-gene case

# You can access this file by (https://www.dropbox.com/scl/fo/b1r27fxqlq84ywory8qmo/ABRuFJpz1FIAy9vluMhLAMU?rlkey=n3pzgkdd0fqz9waamidlyfcfy&dl=0).
# Note!!! The file size is very large.
trans.p = fread('2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt', sep = '\t')

# This matrix includes the p-values from snp -> core genes
load('trans-eQTL_genetype_genename.Rdata')

p.trans = t(pval.sub)

# This file includes the significant cis-gene
load('snp2sig_cisgene.Rdata')

load('GCST90082964.Rdata')
head(dat_M3.1)

out.pheno <- strsplit(as.character(dat_M3.1$Name),'\\.')

info.out.pheno = do.call(rbind, out.pheno)

out.new <- strsplit(as.character(info.out.pheno[,1]),'\\(')

info.out.new = do.call(rbind, out.new)

dat_M3.1$Name = info.out.new[,1]

dat_M3.1_new = dat_M3.1[!dat_M3.1$Name%in%names(which(table(dat_M3.1$Name)!=1)),]


head(dat_M3.1_new)

cau = data.frame(pval = dat_M3.1_new$p_value)
rownames(cau) = dat_M3.1_new$Name


p.wgs <- cau[,1]
p.wgs <- na.omit(p.wgs)
names(p.wgs) = dat_M3.1_new$Name

# Be care about the build version of the position
ref.table = read.delim('gene_position.txt', sep = '\t')

threshold = 0.05/nrow(cau)
eta.wgs = threshold

## only keep gene type in (lincRNA, protein_coding)

ref.table.keep <- ref.table[ref.table$type %in% c('lincRNA', 'protein_coding'),]

## remove genes on chr M, X, and Y
ref.table.keep <- ref.table.keep[!(ref.table.keep$Chromosome %in% c('chrM', 'chrX', 'chrY')),]

ref.table.keep <-ref.table.keep[!duplicated(ref.table.keep$gene_name),]


## standardize gene name
out.ref <- strsplit(as.character(ref.table.keep$gene_name),'\\.')

out.ref.new = do.call(rbind, out.ref)

ref.table.keep$gene_name = out.ref.new[,1]

# Here, we consider all gene1 included in UKBB. You may specify a vector of candidate genes 1

candidate = dput(colnames(p.trans))



# STEP 1: we applied COTA with UKBB and eQTLGen data with target fdr level at 0.1
result = med_gene(p.trans, p.wgs, ref.table, candidate, target.fdr=0.1, dist=5e6, gene1.type = 'SNP', SNP.ref = trans.p)

# STEP 2: we obatined the pair identified by COTA
result.pair = calc_pair.snp(mat.sig = result$mat.sig, mat.p = result$mat.p, p.wgs, ref.table.keep = ref.table.keep, gene1 = colnames(result$mat.sig), 
                            uniq_snp, eta.wgs=1e-5, GRCh = '17')


# STEP 3: we plotted the pair identified by COTA
setwd('/Users/px/Desktop/test/example/')

## load conservation score table
conv.table = read_excel('conservation_scores.xlsx')

## set saving directory of the figures
pic_dir = '/Users/px/Desktop/test/example/'



gen_fig(result.pair$gene.pair,  p.wgs, eta.wgs=1e-5, pic_dir)

```

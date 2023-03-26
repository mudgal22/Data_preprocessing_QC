###Author : Urvy Mudgal 
###Description : RScript for Data preprocessing and quality control 
##Includes Normalization, RLE/NUSE statistics, batch effect correction, and PCA 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)
library(ggfortify)
library(dplyr)

# R runtime
ptm <- proc.time()


file_location = 'path of project file'

celpath <- system.file('celfiles', package='affydata')
fns <- list.celfiles(path=celpath, full.names=TRUE)

# read in files
celbatch <- ReadAffy(celfile.path=file_location)

# normalize files together
rma <- rma(celbatch)

# convert AffyBatch to PLMset
plm_set <- fitPLM(celbatch, normalize = TRUE, background = TRUE)

# relative log expression (RLE)
rle_stats <- data.frame(t(affyPLM::RLE(plm_set, type='stats')))

# plot rle_stats
rle_medians<-ggplot(rle_stats, aes(x=median)) + 
  geom_histogram(bins=50, 
                 color = 'black', 
                 fill = 'maroon') +
  labs(title = 'Medians_RLE') +
  theme_classic()

# normalized unscaled standard errors (NUSE)
nuse_stats <- data.frame(t(NUSE(plm_set, type = 'stats')))

# plot nuse_stats
nuse_medians<- ggplot(nuse_stats, aes(x=median)) + 
  geom_histogram(bins=50, 
                 color = 'black', 
                 fill = 'maroon') +
  labs(title = 'Medians_NUSE') +
  theme_classic()

# correction for batch effects

annotation_csv <- read.csv(file = 'path of metadata.csv file')

edata <- exprs(rma)
batch <- annotation_csv$normalizationcombatbatch

modcombat <- model.matrix(~as.factor(normalizationcombatmod), data = annotation_csv)

combat_edata = ComBat(dat = edata, batch = batch, mod = modcombat)

# write combat_edata to csv
write.csv(combat_edata, file = 'path to save new edata.csv file')

# transpose before pca
trans_edata <- t(combat_edata)

# scale and center
scaled_edata <- scale(trans_edata, center = TRUE, scale = TRUE)

# retranspose
scaled_edata <- t(scaled_edata)


# perform pca
pca <- prcomp(scaled_edata, center = FALSE, scale = FALSE)

# pull variance explained
explained_var <- pca$sdev^2 / sum(pca$sdev^2)
explained_var[1:5]

pca$rotation = as.data.frame(pca$rotation)

# plot pca
pc1vspc2<-pca$rotation %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(size = 0.5) +
  theme_grey() + 
  labs(x = paste0('PC1: ', round(explained_var[1]*100, 2), '%'),
       y = paste0('PC2: ', round(explained_var[2]*100, 2), '%'),
       title = 'PC1 vs PC2')


# R runtime stop
runtime <- proc.time() - ptm
runtime
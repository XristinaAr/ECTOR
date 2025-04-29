## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Autor: Francisco J. Romero-Campero - fran@us.es

## Instalar chipseeker y paquete de anotación de Arabidopsis thaliana
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")

library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28


## Leer fichero de picos
prr5.peaks <- readPeakFile(peakfile = "prr5_peaks.narrowPeak",header=FALSE)

## Definir la región que se considera promotor entorno al TSS
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)

## Anotación de los picos
prr5.peakAnno <- annotatePeak(peak = prr5.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)

plotAnnoPie(prr5.peakAnno)
plotAnnoBar(prr5.peakAnno)
plotDistToTSS(prr5.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")
upsetplot(prr5.peakAnno)

plotPeakProf2(peak = prr5.peaks, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)


## Convertir la anotación a data frame
prr5.annotation <- as.data.frame(prr5.peakAnno)
head(prr5.annotation)
dim(prr5.annotation)

target.genes <- prr5.annotation$geneId[prr5.annotation$annotation == "Promoter"]
vamos a filtrar aquellos genes cuyos picos se encuentran en la zona de su promotor.
write(x = target.genes,file = "prr5_target_genes.txt")


## Enriquecimiento funcional. 


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(org.At.tair.db)

prr5.enrich.go <- enrichGO(gene = target.genes,
                           OrgDb         = org.At.tair.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable      = FALSE,
                           keyType = "TAIR")

barplot(prr5.enrich.go,showCategory = 20)
dotplot(prr5.enrich.go,showCategory = 20)

library(enrichplot)
emapplot(pairwise_termsim(prr5.enrich.go),showCategory = 20)
cnetplot(prr5.enrich.go,showCategory = 6, max.overlaps = 20)
cnetplot(prr5.enrich.go,showCategory = 6, cex_label_category = 1,
         cex_label_gene = 0.5, max.overlaps = Inf )



prr5.enrich.kegg <- enrichKEGG(gene  = target.genes,
                               organism = "ath",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05)
df.prr5.enrich.kegg <- as.data.frame(prr5.enrich.kegg)
head(df.prr5.enrich.kegg)

## ChIPpeakAnno es un paquete de R de Bioconductor que implementa 
## análisis de los resultados del procesamiento de los datos de ChIP-seq
## alternativo y complementario al presentado anteriormente. 

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")

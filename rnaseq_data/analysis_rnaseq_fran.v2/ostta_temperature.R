##########################################################################
## RNA-seq data analysis for temperature response in Ostreococcus tauri ##
##                                                                      ##    
## Authors: Christina Arvanitidou carvanitidou@us.es                    ##
##          Francisco J. Romero-Campero fran@us.es                      ##
##########################################################################

## The ballgown package is used to load the data processed using 
## the mapper hisat2 and the transcripts assembler stringtie
library(ballgown)

## In order to load the data it is necessary to specify the experimental
## design.
experimental.design <- read.csv("experimental_design.csv",stringsAsFactors = T, sep = "\t")
experimental.design

## Load the gene expression estimation data
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=experimental.design)
bg.data
sampleNames(bg.data)

## Extract gene expression data
gene.expression <- gexpr(bg.data)
head(gene.expression)
dim(gene.expression)
gene.names <- rownames(gene.expression)
write(x = gene.names,file = "ostta_temp_gene_names.txt")

## Naming columns according to the experimental design 
colnames(gene.expression) <- c(paste("T10",1:3,sep="_"),
                               paste("T14",1:3,sep="_"),
                               paste("T20",1:3,sep="_"),
                               paste("T26",1:3,sep="_"),
                               paste("T24",1:3,sep="_"))

## Adding 1 to the gene expression matrix so we can apply log2 transformation
## during normalization (log2(0) is not defined)
gene.expression.1 <- gene.expression + 1


## Saving gene expression data
write.table(x = gene.expression.1,file = "temperature_gene_expression.tsv",
            quote = F,row.names = F,
            sep = "\t")


## Boxplot representing global gene expression over all samples
boxplot(gene.expression, outline=F,col=rainbow(12),ylab="Gene Expression (FPKM) + 1",
        cex.lab=1.5,las=2)

boxplot(log2(gene.expression.1), outline=F,col=rainbow(12),
        ylab="log2(Gene Expression + 1)",
        cex.lab=1.5,las=2)

## Data Normalization
library(NormalyzerDE)

## The experimental desing needs to be specified in a different way for
## the package normalyzerDE
design <- data.frame(sample=colnames(gene.expression),
                     group=experimental.design$temperature)

write.table(x = design,file = "normalyzer_design.tsv",quote = F,row.names = F,
            sep = "\t")

## Apply normalization and generate reports with their evaluation
normalyzer(jobName = "temperature_normalization",designPath = "normalyzer_design.tsv",
           dataPath = "temperature_gene_expression.tsv",outputDir = ".", skipAnalysis = T)

## Quantile normalization is chosen as the best method for our data
normalized.gene.expression <- read.table(file="temperature_normalization/Quantile-normalized.txt", header=T)
head(normalized.gene.expression)
rownames(normalized.gene.expression) <- gene.names
 
## Boxplot representing normalized global gene expression over all samples
boxplot(normalized.gene.expression, outline=F,col=rainbow(12),
        ylab="log2(FPKM + 1)",cex.lab=1.5,las=2)


## Principal Component Analysis
library(FactoMineR)
library(factoextra)

## Tutorials:
## http://www.sthda.com/english/articles/tag/factominer/
## http://www.sthda.com/english/
## http://www.sthda.com/english/articles/22-principal-component-methods-videos/65-pca-in-r-using-factominer-quick-scripts-and-videos/
## http://www.sthda.com/english/articles/22-principal-component-methods-videos/74-hcpc-using-factominer-video/

## Gene expression matrix needs to be formatted in a specific way for factominer
pca.gene.expression <- data.frame(colnames(normalized.gene.expression),
                                  t(normalized.gene.expression))
colnames(pca.gene.expression)[1] <- "Sample"
#head(pca.gene.expression)

## Computing PCA and Hierarchical clustering
res.pca <- PCA(pca.gene.expression, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 4)   

## Visualizing HCPC
png(filename = "hierarchical_clustering.png")
fviz_dend(res.hcpc,k=3,main = "Samples Hierarchical Clustering",
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 500      # Augment the room for labels
)
dev.off()

png(filename = "pca.png")
fviz_pca_ind(res.pca, col.ind = experimental.design[,2], 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = T,
             legend.title="Temperature",
             title="",
             show_legend=TRUE,show_guide=TRUE) + scale_color_manual(values=c("#00FFFF","#6699FF", "#CC33FF","#FF0000"))
dev.off()

expression.matrix <- 2^normalized.gene.expression - 1


gene.expression.barplot <- function(gene,expression.matrix)
{
 cond.names <- c("T10ºC","T14ºC","T20ºC","T26ºC")
 gene.expr <- matrix(nrow=length(cond.names),ncol=3)
 rownames(gene.expr) <- cond.names
 for(i in 1:length(cond.names))
 {
  gene.expr[i,] <- unlist(c(expression.matrix[gene, paste(cond.names[i],1:3,sep="_") ]))
 }
 
 means <- apply(X = gene.expr,MARGIN = 1,FUN = mean)
 sds <- apply(X = gene.expr,MARGIN = 1,FUN = sd)
 sds <- sds / 3
 
 arrow.top <- means + sds
 arrow.bottom <- means - sds
 
 
 
 
 xpos <- barplot(means,ylim=c(0,1.2*max(arrow.top)),col=c("#00FFFF","#6699FF", "#CC33FF","#FF0000"), #rainbow(length(cond.names)),
                 main=gene,
                 names.arg = c("10ºC","14ºC","20ºC","26ºC"),cex.main=2,
                 ylab="FPKM",las=2)
 arrows(xpos, arrow.top, xpos, arrow.bottom,code = 3,angle=90,length=0.05)
 box()
}

i <- 1



for(i in 1:length(target.genes))
{
 current.gene <- target.genes[i]
 png(filename = paste0(paste0("gene_barplots/",current.gene),".png"))
 gene.expression.barplot(gene=current.gene, expression.matrix)
 dev.off()
}



gene.expression.barplot(gene="ostta02g01020",expression.matrix)
gene.expression.barplot(gene="ostta02g00380",expression.matrix)

gene.expression.barplot(gene="ostta10g00010", gene.name="ostta10g00010",expression.matrix)

gene.expression.barplot(gene="ostta10g02060", gene.name="ostta10g02060",expression.matrix)


gene.expression.barplot(gene="ostta12g02370", gene.name="ostta02g05010",expression.matrix)



gene.expression.barplot(gene="ostta08g00390", gene.name="?",expression.matrix,cond.names=c("T14","T20","T26"))



## The package limma is used for differential gene expression.
library(limma)

## Limma needs again a different way of specifying the experimental design
limma.experimental.design <- model.matrix(~ -1+factor(c(1,2,3,4,5,6,7,8,9,10,11,12,
                                                        1,2,3,4,5,6,7,8,9,10,11,12,
                                                        1,2,3,4,5,6,7,8,9,10,11,12)))
colnames(limma.experimental.design) <- experimental.design$treatment[1:12]

##Linear model fit (compute average gene expression among replicates)
linear.fit <- lmFit(normalized.gene.expression, limma.experimental.design)

## Contrasts Specification experimental_condition-control_condition
contrast.matrix <- makeContrasts(S14_A-A24_A,
                                 S14_B-A24_B,
                                 S14_C-A24_C,
                                 levels = apply(X = experimental.design[,2:3],
                                                MARGIN = 1,
                                                FUN = function(x)
                                                 {return(paste(x,collapse="_"))})[1:12])

## Fold change and q-value computation
contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## DEGs for chilling 30min.
chilling.30min <- topTable(contrast.results, 
                           number=nrow(normalized.gene.expression),
                           coef=1,sort.by="logFC")
head(chilling.30min)

log.fold.change.chilling.30min <- chilling.30min$logFC
q.value.chilling.30min <- chilling.30min$adj.P.Val
genes.ids.chilling.30min <- rownames(chilling.30min)
names(log.fold.change.chilling.30min) <- genes.ids.chilling.30min
names(q.value.chilling.30min) <- genes.ids.chilling.30min

activated.genes.chilling.30min <- genes.ids.chilling.30min[log.fold.change.chilling.30min > 1 & q.value.chilling.30min < 0.05]
repressed.genes.chilling.30min <- genes.ids.chilling.30min[log.fold.change.chilling.30min < - 1 & q.value.chilling.30min < 0.05]

length(activated.genes.chilling.30min)
length(repressed.genes.chilling.30min)

write(x = activated.genes.chilling.30min,file = "activated_S14_vs_A24_30min.txt")
write(x = repressed.genes.chilling.30min,file = "repressed_S14_vs_A24_30min.txt")


## Volcano plot
log.q.val.chilling.30min <- -log10(q.value.chilling.30min)
plot(log.fold.change.chilling.30min,
     log.q.val.chilling.30min,
     pch=19,col="grey",cex=0.8,
     xlim=c(-4,4),ylim = c(0,4), 
     xlab="log2(Fold-chage)",ylab="-log10(q-value)",cex.lab=1.5)

points(x = log.fold.change.chilling.30min[activated.genes.chilling.30min],
       y = log.q.val.chilling.30min[activated.genes.chilling.30min],col="red",cex=0.8,pch=19)
points(x = log.fold.change.chilling.30min[repressed.genes.chilling.30min],
       y = log.q.val.chilling.30min[repressed.genes.chilling.30min],col="blue",cex=0.8,pch=19)


## DEGs for chilling 4h.
chilling.4h <- topTable(contrast.results, 
                           number=nrow(normalized.gene.expression),
                           coef=2,sort.by="logFC")
head(chilling.4h)

log.fold.change.chilling.4h <- chilling.4h$logFC
q.value.chilling.4h <- chilling.4h$adj.P.Val
genes.ids.chilling.4h <- rownames(chilling.4h)
names(log.fold.change.chilling.4h) <- genes.ids.chilling.4h
names(q.value.chilling.4h) <- genes.ids.chilling.4h

activated.genes.chilling.4h <- genes.ids.chilling.4h[log.fold.change.chilling.4h > 1 & q.value.chilling.4h < 0.05]
repressed.genes.chilling.4h <- genes.ids.chilling.4h[log.fold.change.chilling.4h < - 1 & q.value.chilling.4h < 0.05]

length(activated.genes.chilling.4h)
length(repressed.genes.chilling.4h)

write(x = activated.genes.chilling.4h,file = "activated_S14_vs_A24_4h.txt")
write(x = repressed.genes.chilling.4h,file = "repressed_S14_vs_A24_4h.txt")


## Volcano plot
log.q.val.chilling.4h <- -log10(q.value.chilling.4h)
plot(log.fold.change.chilling.4h,
     log.q.val.chilling.4h,
     pch=19,col="grey",cex=0.8,
     xlim=c(-6,6),ylim = c(0,10), 
     xlab="log2(Fold-chage)",ylab="-log10(q-value)",cex.lab=1.5)

points(x = log.fold.change.chilling.4h[activated.genes.chilling.4h],
       y = log.q.val.chilling.4h[activated.genes.chilling.4h],col="red",cex=0.8,pch=19)
points(x = log.fold.change.chilling.4h[repressed.genes.chilling.4h],
       y = log.q.val.chilling.4h[repressed.genes.chilling.4h],col="blue",cex=0.8,pch=19)


## DEGs for chilling 24h.
chilling.24h <- topTable(contrast.results, 
                         number=nrow(normalized.gene.expression),
                         coef=3,sort.by="logFC")
head(chilling.24h)

log.fold.change.chilling.24h <- chilling.24h$logFC
q.value.chilling.24h <- chilling.24h$adj.P.Val
genes.ids.chilling.24h <- rownames(chilling.24h)
names(log.fold.change.chilling.24h) <- genes.ids.chilling.24h
names(q.value.chilling.24h) <- genes.ids.chilling.24h

activated.genes.chilling.24h <- genes.ids.chilling.24h[log.fold.change.chilling.24h > 1 & q.value.chilling.24h < 0.05]
repressed.genes.chilling.24h <- genes.ids.chilling.24h[log.fold.change.chilling.24h < - 1 & q.value.chilling.24h < 0.05]

length(activated.genes.chilling.24h)
length(repressed.genes.chilling.24h)

write(x = activated.genes.chilling.24h,file = "activated_S14_vs_A24_24h.txt")
write(x = repressed.genes.chilling.24h,file = "repressed_S14_vs_A24_24h.txt")


## Volcano plot
log.q.val.chilling.24h <- -log10(q.value.chilling.24h)
plot(log.fold.change.chilling.24h,
     log.q.val.chilling.24h,
     pch=19,col="grey",cex=0.8,
     xlim=c(-6,6),ylim = c(0,8), 
     xlab="log2(Fold-chage)",ylab="-log10(q-value)",cex.lab=1.5)

points(x = log.fold.change.chilling.24h[activated.genes.chilling.24h],
       y = log.q.val.chilling.24h[activated.genes.chilling.24h],col="red",cex=0.8,pch=19)
points(x = log.fold.change.chilling.24h[repressed.genes.chilling.24h],
       y = log.q.val.chilling.24h[repressed.genes.chilling.24h],col="blue",cex=0.8,pch=19)
























## Enriquecimiento funcional. 
library(clusterProfiler)
library(org.At.tair.db)

activated.atha.enrich.go <- enrichGO(gene          = activated.genes,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

barplot(activated.atha.enrich.go,showCategory = 20)
dotplot(activated.atha.enrich.go,showCategory = 20)
emapplot(activated.atha.enrich.go,showCategory = 20)
cnetplot(activated.atha.enrich.go,showCategory = 20)


repressed.atha.enrich.go <- enrichGO(gene          = repressed.genes,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

barplot(repressed.atha.enrich.go,showCategory = 20)
dotplot(repressed.atha.enrich.go,showCategory = 20)
emapplot(repressed.atha.enrich.go,showCategory = 20)
cnetplot(repressed.atha.enrich.go,showCategory = 20)


activated.atha.enrich.kegg <- enrichKEGG(gene  = activated.genes,
                                         organism = "ath",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05)
df.activated.atha.enrich.kegg <- as.data.frame(activated.atha.enrich.kegg)
head(df.activated.atha.enrich.kegg)


repressed.atha.enrich.kegg <- enrichKEGG(gene  = repressed.genes,
                                         organism = "ath",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05)

df.repressed.atha.enrich.kegg <- as.data.frame(repressed.atha.enrich.kegg)
head(df.repressed.atha.enrich.kegg)

## Podemos visualizar las rutas KEGG usando el paquete pathveiw
library(pathview)
pathview(gene.data = sort(log.fold.change,decreasing = TRUE),
         pathway.id = "ath00941",
         species = "ath",
         limit = list(gene=max(abs(log.fold.change)), cpd=1),gene.idtype = "TAIR")




## Código para desarrollar una función gráfico de barras

gene <- "ostta06g02940"
gene.name <- "ostta06g02940"
gene <- "ostta03g04010"

colnames(expression.matrix)




cond.names <- c("A24_B","S14_B","S14_A", "A14_A")
cond.names <- apply(X = experimental.design[,2:3],
                        MARGIN = 1,
                        FUN = function(x){return(paste(x,collapse="_"))})[1:12]


gene.name <- beta

gene.expression.barplot <- function(gene, gene.name,expression.matrix,cond.names)
{
 gene.expr <- matrix(nrow=length(cond.names),ncol=3)
 for(i in 1:length(cond.names))
 {
  gene.expr[i,] <- unlist(c(expression.matrix[gene, paste(cond.names[i],1:3,sep="_") ]))
 }
 
 means <- apply(X = gene.expr,MARGIN = 1,FUN = mean)
 sds <- apply(X = gene.expr,MARGIN = 1,FUN = sd)
 sds <- sds / 3

  arrow.top <- means + sds
  arrow.bottom <- means - sds
  
  
  

  xpos <- barplot(means,ylim=c(0,1.2*max(arrow.top)),col=c(rep("blue",3),rep("darkblue",3),rep("red",3),rep("darkred",3)), #rainbow(length(cond.names)),
                  main=gene.name,
                  names.arg = cond.names,cex.main=2,
                  ylab="FPKM",las=2)
  arrows(xpos, arrow.top, xpos, arrow.bottom,code = 3,angle=90,length=0.05)
}
text.main <- bquote(beta - .(gene))

all.conds <- apply(X = experimental.design[,2:3],
                   MARGIN = 1,
                   FUN = function(x){return(paste(x,collapse="_"))})[1:12]
all.conds <- c("S14_A", "S14_B", "S14_C", "A14_A", "A14_B", "A14_C",
               "S24_A", "S24_B", "S24_C", "A24_A", "A24_B", "A24_C")


expression.matrix <- 2^normalized.gene.expression - 1
gene.expression.barplot(gene="ostta06g02940", gene.name="ostta06g02940",expression.matrix,cond.names=all.conds)

gene.expression.barplot(gene="ostta06g02340", expression.matrix,cond.names=all.conds)


tfs <- read.table(file = "transcription_factors_list.tsv",sep="\t",header=T)
head(tfs)

i <- 1
for(i in 1:nrow(tfs))
{
 tf.name <- paste(tfs[i,1],tfs[i,2],sep="_")
 jpeg(filename = paste(tf.name,"jpg",sep="."))
 gene.expression.barplot(gene=tfs[i,2],gene.name = tf.name, 
                         expression.matrix = expression.matrix,
                         cond.names=all.conds)
 dev.off()
}


 gene.expression.barplot(gene="ostta03g04010",
                        gene.name = "FabZ", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)


gene.expression.barplot(gene="ostta01g03280",
                        gene.name = "ostta01g03280", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)


gene.expression.barplot(gene="ostta12g01000",
                        gene.name = "ostta12g01000", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)


gene.expression.barplot(gene="ostta15g01470",
                        gene.name = "ostta15g01470", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)

gene.expression.barplot(gene="ostta03g01370",
                        gene.name = "ostta03g01370", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)

gene.expression.barplot(gene="ostta02g01840",
                        gene.name = "ostta02g01840", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)



gene.expression.barplot(gene="ostta06g02930",
                        gene.name = "PGM", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)


gene.expression.barplot(gene="ostta07g03440",
                        gene.name = "APL", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)

gene.expression.barplot(gene="ostta06g02940",
                        gene.name = "GBSS", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)


gene.expression.lineplot <- function(gene, gene.name,expression.matrix)
{
 cond.names <- c("S14_A", "S14_B", "S14_C", "A14_A", "A14_B", "A14_C",
                 "S24_A", "S24_B", "S24_C", "A24_A", "A24_B", "A24_C")
 
 gene.expr <- matrix(nrow=length(cond.names),ncol=3)
 for(i in 1:length(cond.names))
 {
  gene.expr[i,] <- unlist(c(expression.matrix[gene, paste(cond.names[i],1:3,sep="_") ]))
 }
 
 means <- apply(X = gene.expr,MARGIN = 1,FUN = mean)
 sds <- apply(X = gene.expr,MARGIN = 1,FUN = sd)
 
 arrow.top <- means + sds/sqrt(3)
 arrow.bottom <- means - sds/sqrt(3)
 
 upper <- ceiling(1.1*max(arrow.top))
 lower <- floor(0.9*min(arrow.bottom))
 
 photoperiod.lane <- floor((upper - lower) / 50)
 
 plot(x=c(0,4,24),means[1:3],type="o",lwd=4,col="blue",axes=F,xlab="",ylab="FPKM",cex.lab=1.5,ylim=c(lower - 4 * photoperiod.lane,upper),main=gene.name,cex.main=2)
 arrows(c(0,4,24), arrow.top[1:3], c(0,4,24), arrow.bottom[1:3],code = 3,angle=90,length=0.05,lwd=4,col="blue")
 
 lines(x=c(0,4,24),means[4:6],type="o",lwd=4,col="darkblue")
 arrows(c(0,4,24), arrow.top[4:6], c(0,4,24), arrow.bottom[4:6],code = 3,angle=90,length=0.05,lwd=4,col="darkblue")
 
 lines(x=c(0,4,24),means[7:9],type="o",lwd=4,col="red")
 arrows(c(0,4,24), arrow.top[7:9], c(0,4,24), arrow.bottom[7:9],code = 3,angle=90,length=0.05,lwd=4,col="red")
 
 lines(x=c(0,4,24),means[10:12],type="o",lwd=4,col="darkred")
 arrows(c(0,4,24), arrow.top[10:12], c(0,4,24), arrow.bottom[10:12],code = 3,angle=90,length=0.05,lwd=4,col="darkred")


 polygon(x = c(0,16,16,0),y=c(lower - 4 * photoperiod.lane, lower - 4 * photoperiod.lane,
                              lower - 0.5 * photoperiod.lane, lower - 0.5 * photoperiod.lane),lwd=4,border="blue")
 polygon(x = c(16,24,24,16),y=c(lower - 4 * photoperiod.lane, lower - 4 * photoperiod.lane,
                                lower - 0.5 * photoperiod.lane, lower - 0.5 * photoperiod.lane),lwd=4,border="blue",col="blue")
 
  
 axis(side=2,lwd=3,las=2,cex.axis=1.5) 
}

for(i in 1:nrow(tfs))
{
 tf.name <- paste(tfs[i,1],tfs[i,2],sep="_")
 jpeg(filename = paste(tf.name,"jpg",sep="."))
 gene.expression.lineplot(gene=tfs[i,2],gene.name = tf.name, 
                         expression.matrix = expression.matrix)
 dev.off()
}


gene.expression.lineplot(gene="ostta11g02780",
                         gene.name = "PGI", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta06g02930",
                        gene.name = "PGM", 
                        expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta07g03440",
                        gene.name = "APL", 
                        expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta07g03070",
                         gene.name = "APS", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta06g02940",
                        gene.name = "GBSS", 
                        expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta03g00870",
                         gene.name = "SBE", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta01g03280",
                         gene.name = "acetyl-CoA carboxylase", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta10g03450",
                         gene.name = "acetyl-CoA carboxylase 2", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta10g01690",
                         gene.name = "malonyl_CoA_methylmalonyl_CoA_synthetase", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta10g01690",
                         gene.name = "malonyl_CoA_methylmalonyl_CoA_synthetase", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta09g03060",
                         gene.name = "FabD1", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta20g00030",
                         gene.name = "FabD2", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta15g01890",
                         gene.name = "FabH", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta02g01840",
                         gene.name = "FabF1", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta02g04460",
                         gene.name = "FabF2", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta15g01470",
                         gene.name = "FabF3", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta01g02220",
                         gene.name = "FabG", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta03g04010",
                         gene.name = "FabZ", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta07g03880",
                         gene.name = "FabI", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta03g01370",
                         gene.name = "MECR1", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta17g00300",
                         gene.name = "MECR2", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta04g03510",
                         gene.name = "ostta04g03510", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta01g06660",
                         gene.name = "ostta01g06660", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta05g00250",
                         gene.name = "ostta05g00250", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta12g01000",
                         gene.name = "ostta12g01000", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta16g00040",
                         gene.name = "ostta16g00040", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta03g01370",
                         gene.name = "ostta03g01370", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta17g00300",
                         gene.name = "ostta17g00300", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta17g00650",
                         gene.name = "ostta17g00650", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta18g00810",
                         gene.name = "ostta18g00810", 
                         expression.matrix = expression.matrix)



gene.expression.lineplot(gene="ostta17g01290",
                         gene.name = "ostta17g01290", 
                         expression.matrix = expression.matrix)


gene.expression.lineplot(gene="ostta02g03030",
                         gene.name = "CSP", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta02g00450",
                         gene.name = "ostta02g00450", 
                         expression.matrix = expression.matrix)

gene.expression.lineplot(gene="ostta04g03630",
                         gene.name = "ostta04g03630", 
                         expression.matrix = expression.matrix)



gene.expression.lineplot(gene="ostta03g03040",
                         gene.name = "ostta03g03040", 
                         expression.matrix = expression.matrix)







gene.expression.lineplot(gene="ostta09g03210",
                         gene.name = "FAD4", 
                         expression.matrix = expression.matrix)


## Fatty acid control  Florence Corellou

png(filename = "FAD4_ostta09g03210.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta09g03210",
                        gene.name = "FAD4", 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "MDGG_D9_D6_ostta17g02260.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta17g02260",
                        gene.name = bquote(MDGG - Delta~9  / Delta~6  - ostta17g02260), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "D4p_ostta13g01550.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta13g01550",
                        gene.name = bquote(Delta~4~p - ostta13g01550), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "sad_D9_ostta04g03510.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta04g03510",
                        gene.name = bquote(SAD~Delta~9 - ostta04g03510), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "D12_1_ostta02g05510.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta02g05510",
                        gene.name = bquote(Delta~12~1 - ostta02g05510), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "omega3_D15_ostta03g03040.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta03g03040",
                        gene.name = bquote(omega~3~~Delta~15 - ostta03g03040), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "D5_ostta12g02740.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta12g02740",
                        gene.name = bquote(Delta~5 - ostta12g02740), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "D5_ostta12g02740.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta12g02740",
                        gene.name = bquote(Delta~5 - ostta12g02740), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "D6_D8_ostta05g00100.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta05g00100",
                        gene.name = bquote( Delta~6 / Delta~8 - ostta05g00100), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "D6_D8_ostta10g02580.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta10g02580",
                        gene.name = bquote( Delta~6 / Delta~8 - ostta10g02580), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "D6_D8_ostta15g01140.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta15g01140",
                        gene.name = bquote( Delta~6 / Delta~8 - ostta15g01140), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "D6_D8_ostta18g01930.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta18g01930",
                        gene.name = bquote( Delta~6 / Delta~8 - ostta18g01930), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "Iwane_ostta12g01890.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta12g01890",
                        gene.name = bquote( Iwane - ostta12g01890), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

## Some examples found in our analysis
png(filename = "ACACA_1_ostta01g03280.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta01g03280",
                        gene.name = bquote( ACACA - ostta01g03280), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "ACACA_2_ostta10g03450.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta10g03450",
                        gene.name = bquote( ACACA - ostta10g03450), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "ACSF3_ostta10g01690.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta10g01690",
                        gene.name = bquote( ACSF3 - ostta10g01690), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "FabD_ostta09g03060.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta09g03060",
                        gene.name = bquote( FabD - ostta09g03060), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "FabD_ostta20g00030.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta20g00030",
                        gene.name = bquote( FabD - ostta20g00030), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "FabF_ostta02g01840.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta02g01840",
                        gene.name = bquote( FabF - ostta02g01840), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "FabF_ostta02g04460.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta02g04460",
                        gene.name = bquote( FabF - ostta02g04460), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "FabF_ostta15g01470.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta15g01470",
                        gene.name = bquote( FabF - ostta15g01470), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "FabH_ostta15g01890.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta15g01890",
                        gene.name = bquote( FabH - ostta15g01890), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "FabG_ostta01g02220.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta01g02220",
                        gene.name = bquote( FabG - ostta01g02220), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "FabZ_ostta03g04010.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta03g04010",
                        gene.name = bquote( FabZ - ostta03g04010), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "FabI_ostta07g03880.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta07g03880",
                        gene.name = bquote( FabI - ostta07g03880), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "MECR_ostta03g01370.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta03g01370",
                        gene.name = bquote( MECR - ostta03g01370), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "MECR_ostta17g00300.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta17g00300",
                        gene.name = bquote( MECR - ostta17g00300), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "desA1_ostta04g03510.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta04g03510",
                        gene.name = bquote( desA1 - ostta04g03510), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "FATA_ostta01g06660.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta01g06660",
                        gene.name = bquote( FATA - ostta01g06660), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "ACSL_ostta05g00250.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta05g00250",
                        gene.name = bquote( ACSL - ostta05g00250), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

## Biosynthesis of unsaturated fatty acids
png(filename = "FADS2_ostta13g01040",width = 400,height = 400)
gene.expression.barplot(gene="ostta13g01040",
                        gene.name = bquote( Delta~6 - ostta13g01040), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "HACD_ostta05g03730.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta05g03730",
                        gene.name = bquote( HACD - ostta05g03730), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "HACD_ostta13g02530.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta13g02530",
                        gene.name = bquote( HACD - ostta13g02530), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "D9_ostta04g03510.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta04g03510",
                        gene.name = bquote( Delta~9 - ostta04g03510), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "D9_ostta17g02010.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta17g02010",
                        gene.name = bquote( Delta~9 - ostta17g02010), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "D8_ostta13g01040.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta13g01040",
                        gene.name = bquote( Delta~8 - ostta13g01040), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "ELOVL4_ostta02g00200.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta02g00200",
                        gene.name = bquote( ELOVL4 - ostta02g00200), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "ELOVL4_ostta09g02530.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta09g02530",
                        gene.name = bquote( ELOVL4 - ostta09g02530), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "PPT_ostta17g00650.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta17g00650",
                        gene.name = bquote( PPT - ostta17g00650), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "TER_ostta03g04210.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta03g04210",
                        gene.name = bquote( TER - ostta03g04210), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()


png(filename = "ACSL_ostta12g01000.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta12g01000",
                        gene.name = bquote( ACSL - ostta12g01000), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "ACSL_ostta16g00040.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta16g00040",
                        gene.name = bquote( ACSL - ostta16g00040), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "ACOX_ostta01g02650.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta01g02650",
                        gene.name = bquote( ACOX - ostta01g02650), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "ACOX_ostta12g02310.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta12g02310",
                        gene.name = bquote( ACOX - ostta12g02310), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "ACADM_ostta01g00320.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta01g00320",
                        gene.name = bquote( ACADM - ostta01g00320), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "MFP2_ostta03g04620.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta03g04620",
                        gene.name = bquote( MFP2 - ostta03g04620), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "MFP2_ostta09g03720.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta09g03720",
                        gene.name = bquote( MFP2 - ostta09g03720), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "ACAA1_ostta06g04450.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta06g04450",
                        gene.name = bquote( ACAA1 - ostta06g04450), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "SGPL1_ostta02g02240.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta02g02240",
                        gene.name = bquote( SGPL1 - ostta02g02240), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "SPT_ostta16g01250.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta16g01250",
                        gene.name = bquote( SPT - ostta16g01250), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "lacZ_ostta07g03640.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta07g03640",
                        gene.name = bquote( lacZ - ostta07g03640), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()

png(filename = "SMPD2_ostta07g00630.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta07g00630",
                        gene.name = bquote( SMPD2 - ostta07g00630), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()



png(filename = "SMPD2_ostta07g00630.png",width = 400,height = 400)
gene.expression.barplot(gene="ostta02g05010",
                        gene.name = bquote( SMPD2 - ostta07g00630), 
                        expression.matrix = expression.matrix,
                        cond.names=all.conds)
dev.off()




## Análisis de expresión génica diferencial con DESeq2

library(DESeq2)

pheno.data

gene.count.matrix <- read.table(file = "gene_count_matrix.csv",header = T,sep = ",")
head(gene.count.matrix)


gene.ids <- gene.count.matrix$gene_id

gene.count.matrix <- gene.count.matrix[,-1]
rownames(gene.count.matrix) <- gene.ids
head(gene.count.matrix)
head(experimental.design)
nrow(experimental.design)
ncol(gene.count.matrix)
dds <- DESeqDataSetFromMatrix(countData=gene.count.matrix, colData=experimental.design, design = ~ treatment)




#dds$treatment <- relevel(dds$treatment, ref = "A24_B")

dds <- DESeq(dds) 
resultsNames(dds)

## S14_A vs A24_A  
res <- results(dds, contrast = c("treatment","S14_A","A24_A")) 
res

log.fold.change <- res$log2FoldChange
q.value <- res$padj
names(log.fold.change) <- rownames(res)
names(q.value) <- rownames(res)

activated.genes.deseq2 <- rownames(res)[log.fold.change > 1 & q.value < 0.05]
activated.genes.deseq2 <- activated.genes.deseq2[!is.na(activated.genes.deseq2)]
length(activated.genes.deseq2)

repressed.genes.deseq2 <- rownames(res)[log.fold.change < - 1 & q.value < 0.05]
repressed.genes.deseq2 <- repressed.genes.deseq2[!is.na(repressed.genes.deseq2)]

length(activated.genes.deseq2)
length(repressed.genes.deseq2)

write(x = activated.genes.deseq2,file = "activated_S14_A_vs_A24_A.txt")
write(x = repressed.genes.deseq2,file = "repressed_S14_A_vs_A24_A.txt")

tfs.df <- read.table(file = "transcription_factors_list.tsv",sep="\t",header=T)
head(tfs.df)
tfs <- tfs.df$gene_id

tfs.activated <- intersect(activated.genes.deseq2,tfs)
tfs.repressed <- intersect(repressed.genes.deseq2,tfs)

subset(tfs.df, gene_id %in% tfs.activated)
subset(tfs.df, gene_id %in% tfs.repressed)

activated.atha.enrich.go <- enrichGO(gene          = activated.genes.deseq2,
                                     OrgDb         = org.At.tair.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "TAIR")

barplot(activated.atha.enrich.go,showCategory = 20)
dotplot(activated.atha.enrich.go,showCategory = 20)
emapplot(activated.atha.enrich.go,showCategory = 20)
cnetplot(activated.atha.enrich.go,showCategory = 20)


library("TxDb.Otauri.JGI")

install.packages("/home/fran/Dropbox/github_repos/AlgaeFUN/packages/txdb_packages/TxDb.Otauri.JGI/",repos=NULL,type="source")
txdb <- TxDb.Otauri.JGI

my.genes <- c("ostta05g00100",
              "ostta12g02740",
              "ostta15g01140",
              "ostta04g03510",
              "ostta03g03040",
              "ostta10g01690",
              "ostta13g01040"
)

genes.txdb <- as.data.frame()

as.data.frame(resize(subset(genes(txdb), gene_id %in% my.genes),width=1,fix="start"))[,1:3]

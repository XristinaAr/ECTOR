library(seqinr)

ostta.genome <- read.fasta(file = "ostta_genome/ostreococcus_tauri.fa",seqtype = "DNA")
chr.lens <- sapply(X = getSequence(ostta.genome),FUN = length)

number.chrs <- length(getSequence(ostta.genome))
max.chr.len <- max(chr.lens)
chr.lens <- rev(chr.lens)

jpeg(filename = "images/h3k27me3_genome_wide_distribution.jpg",res = 300,width = 900,height = 900)
par(mar=c(0,0,0,0))
plot(x=c(0,max.chr.len), y=c(0,4*number.chrs), col = "white", xlab = "", ylab = "", axes=F)  

for(i in 1:number.chrs)
{
 polygon(x = c(0, chr.lens[i], chr.lens[i], 0),   
         y = c(4*i+1, 4*i+1, 4*i-1, 4*i-1),    
         col = "white") 
}

ld.zt8.1 <- read.table(file="chipseq_data/bed_data/h3k27me3/temp20c/h3k27me3_ld_20_zt8_chip_1_counts.bed",header = F)
ll.20 <- read.table(file="chipseq_data/bed_data/h3k27me3/h3k27me3_ll_20.bed",header = F)
ll.20 <- ll.20[,1:3]
ld.zt8.2 <- read.table(file="chipseq_data/bed_data/h3k27me3/temp20c/h3k27me3_ld_20_zt8_chip_2_counts.bed",header = F)
nrow(ld.zt8.1)

jpeg(filename = "images/h3k27me3_genome_wide_distribution_ll_20.jpg",res = 300,width = 900,height = 900)
par(mar=c(0,0,0,0))
plot(x=c(0,max.chr.len), y=c(0,4*number.chrs), col = "white", xlab = "", ylab = "", axes=F)  

for(i in 1:number.chrs)
{
 polygon(x = c(0, chr.lens[i], chr.lens[i], 0),   
         y = c(4*i+1, 4*i+1, 4*i-1, 4*i-1),    
         col = "white") 
}

for(i in 1:nrow(ll.20))
{
 current.chr <- ll.20[i,1]
 current.start <- ll.20[i,2]
 current.end <- ll.20[i,3]
 
 current.line <- (20:1)[current.chr]
 
 polygon(x = c(current.start, current.end, current.end, current.start),   
         y = c(4*current.line+1, 4*current.line+1, 4*current.line-1, 4*current.line-1),    
         col = "blue",border="blue")
}
dev.off()


for(i in 1:nrow(ll.20))
{
 current.chr <- ld.zt8.[i,1]
 current.start <- ld.zt8.1[i,2]
 current.end <- ld.zt8.1[i,3]
# current.count <- (ld.zt8.1[i,4] + ld.zt8.2[i,4])/2
 
 current.line <- (20:1)[current.chr]
 
 polygon(x = c(current.start, current.end, current.end, current.start),   
         y = c(4*current.line+1, 4*current.line+1, 4*current.line-1, 4*current.line-1),    
         col = "blue",border="blue")
}
dev.off()


h3k27me3.peaks.len <- ll.20$V3 - ll.20$V2
sum(h3k27me3.peaks.len)/sum(chr.lens)
ll.20 <- cbind(ll.20,h3k27me3.peaks.len)
colnames(ll.20) <- c("chr","start","end","length")
head(ll.20)
number.of.h3k27.peaks.chr <- vector(mode = "numeric",length = 20)
peaks.len.chr <- vector(mode = "numeric",length = 20)

for(i in 1:20)
{
 number.of.h3k27.peaks.chr[i] <- nrow(subset(ll.20, chr == i))
 peaks.len.chr[i] <- sum(subset(ll.20, chr == i)[,4]) 
}

percentage.chr.k27 <- peaks.len.chr/rev(chr.lens)

z.scores <- (percentage.chr.k27 - mean(percentage.chr.k27))/sd(percentage.chr.k27)

which(z.scores > 2.5)

p.vals <- pnorm((percentage.chr.k27 - mean(percentage.chr.k27[-19]))/sd(percentage.chr.k27[-19]),lower.tail = F)

values.for.colors <- round(-log10(p.vals)*40)
values.for.colors[19] <- 100

white.blue.palette <- colorRampPalette(c('white','blue'))

jpeg(filename = "images/h3k27me3_chr_barplot_ll_20.jpg",res = 300,width = 1800,height = 900)
barplot(100*percentage.chr.k27,border="black",names.arg = paste0("Chr",1:20),las=2,col=white.blue.palette(100)[values.for.colors]
)
abline(h = 100*mean(percentage.chr.k27[-19]),lty=5)
dev.off()


## Gene features
library(ChIPseeker)
library(TxDb.Otauri.JGI)
txdb <- TxDb.Otauri.JGI
library(org.Otauriv4.eg.db)

ostta.genes <- as.data.frame(genes(txdb))
head(ostta.genes)
i <- 1
i <- 19
#7 y 8 intergenic

target.genes <- vector(mode = "character",length = nrow(ld.zt8.1))
peak.annotation <- vector(mode = "character",length = nrow(ld.zt8.1))

for(i in 1:nrow(ld.zt8.1))
{
 print(i)
 current.chr <- ld.zt8.1[i,1]
 current.start <- ld.zt8.1[i,2]
 current.end <- ld.zt8.1[i,3]

 complete.gene.body <- subset(ostta.genes, seqnames == current.chr & start >= current.start & end <= current.end)
 overlap.tss.or.tes <- subset(ostta.genes, seqnames == current.chr & 
                               start >= current.start & start <= current.end & end > current.end)
 overlap.tes.or.tss <- subset(ostta.genes, seqnames == current.chr & 
                               start < current.start & end >= current.start & end <= current.end)
 inside.gene.body <- subset(ostta.genes, seqnames == current.chr & 
                             start <= current.start & end >= current.end)
 
 
 if(nrow(complete.gene.body) > 0)
 {
  target.genes[i] <- paste(unlist(complete.gene.body$gene_id),collapse = ",")
  peak.annotation[i] <- "gene.body"
 } else if(nrow(overlap.tss.or.tes) > 0)
 {
  if(overlap.tss.or.tes$strand == "+")
  {
   target.genes[i] <- unlist(overlap.tss.or.tes$gene_id)
   peak.annotation[i] <- "TSS"
  } else if (nrow(overlap.tes.or.tss) > 0)
  {
   if(overlap.tes.or.tss$strand == "-")
   {
    target.genes[i] <- unlist(overlap.tes.or.tss$gene_id)
    peak.annotation[i] <- "TSS"
   } else
   {
    target.genes[i] <- unlist(overlap.tss.or.tes$gene_id)
    peak.annotation[i] <- "TES"
   }
  } else
  {
   target.genes[i] <- unlist(overlap.tss.or.tes$gene_id)
   peak.annotation[i] <- "TES"
  }
 } else if(nrow(overlap.tes.or.tss) > 0)
 {
  if(overlap.tes.or.tss$strand == "+")
  {
   target.genes[i] <- unlist(overlap.tes.or.tss$gene_id)
   peak.annotation[i] <- "TES"
  } else if (overlap.tes.or.tss$strand == "-")
  {
   target.genes[i] <- unlist(overlap.tes.or.tss$gene_id)
   peak.annotation[i] <- "TSS"
  }
 } else if(nrow(inside.gene.body) > 0)
 {
  target.genes[i] <- unlist(inside.gene.body$gene_id)
  peak.annotation[i] <- "inner.gene.body"
 } else
 {
  target.genes[i] <- ""
  peak.annotation[i] <- "intergenic"
 }
}

sum(table(peak.annotation))
pie(table(peak.annotation))
as.vector(target.genes)
length(unlist(unique(target.genes)))
is.vector(target.genes)


target.genes[peak.annotation == "inner.gene.body"]


h3k27me3.peaks <- readPeakFile(peakfile = "chipseq_data/bed_data/h3k27me3/temp20c/h3k27me3_ld_20_zt8_chip_1_counts.bed")
covplot(h3k27me3.peaks, weightCol="V4")

promoter <- getPromoters(TxDb=txdb, upstream=200, downstream=200)

peakAnno <- annotatePeak("chipseq_data/bed_data/h3k27me3/temp20c/h3k27me3_ld_20_zt8_chip_1_counts.bed", 
                         tssRegion=c(-200, 200),
                         TxDb=txdb, 
                         annoDb="org.Otauriv4.eg.db")


h3k27me3.annotation <- as.data.frame(peakAnno)

h3k27me3.annotation[h3k27me3.annotation$geneId == "ostta01g01610",]

k27.genes <- unique(subset(h3k27me3.annotation, distanceToTSS == 0)$geneId)
length(k27.genes)
length(k27.genes)/7668

k27.genes <- unique(subset(h3k27me3.annotation, annotation == "Promoter")$geneId)


k27.genes <- unique(subset(h3k27me3.annotation, annotation == "Promoter")$geneId)


h3k27me3.annotation$annotation

write(x = k27.genes,file = "k27_genes.txt")



plotAnnoPie(peakAnno)

jpeg(filename = "images/h3k27me3_pie_chart_ll_20.jpg",res = 300,width = 900,height = 900)
par(mar=c(0,0,0,0))
pie(freq.peak.annotation[c("gene.body",
                           "internal.gene.body",
                           "TSS", 
                           "TES",
                           "intergenic")],labels = c("Full Gene Body", "Gene Body", "TSS", "TES", "Intergenic"),
    clockwise = T,border = "black",col=c("blue","cyan","green","red","grey"))

dev.off()

help(pie)


target.genes <- (unlist(sapply(X = target.genes,FUN = function(x) {strsplit(x,",")[[1]]})))
names(target.genes) <- NULL
target.genes <- unique(target.genes)
write(x = target.genes,file = "target_genes_ll_20.txt")

## metagene plots
library(ChIPpeakAnno)
library(rtracklayer)
library(TxDb.Otauri.JGI)

txdb <- TxDb.Otauri.JGI


k27.gene.body <- (unlist(sapply(X = subset(ll.20.peak.annotation, peak.annotation == "gene.body")[["target.genes"]],
                                FUN = function(x) {strsplit(x,",")[[1]]})))
names(k27.gene.body) <- NULL

## Specification of the set of genes to consider to construct the metageneplot
genes <- k27.gene.body
#genes <- c(h2.k27) ## Example for the set of genes with both marks H2AK121ub and H3K27me3
length(genes)



bw.file <- import.bw(con = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw")
bw.file <- as.data.frame(bw.file)
head(bw.file)

i <- 1


current.gene.chr <- ostta.genes[genes[i],"seqnames"]
current.gene.start <- ostta.genes[genes[i],"start"]
current.gene.end <- ostta.genes[genes[i],"end"]
current.gene.strand <- ostta.genes[genes[i],"strand"]

gene.body.signal <- subset(bw.file, seqnames == current.gene.chr & 
                            start >= current.gene.start & 
                            end <= current.gene.end)[["score"]]

plot(gene.body.signal,type="l")


upstream.start.signal <- subset(bw.file, seqnames == current.gene.chr & 
                                 start >= (current.gene.start - 300) & 
                                 end <= current.gene.start)[["score"]]

plot(upstream.start.signal,type="l")

downstream.end.signal <- subset(bw.file, seqnames == current.gene.chr & 
                                 start >= current.gene.end & 
                                 end <= (current.gene.end + 300))[["score"]]

plot(downstream.end.signal,type="l")

plot(c(upstream.start.signal,gene.body.signal,downstream.end.signal),type="l")

signal.vector <- gene.body.signal

number.tiles <- 50

ntile.signal <- function(signal.vector,number.tiles)
{
 step.len <- floor(length(signal.vector)/number.tiles)
 number.plus.1 <- length(signal.vector) %% number.tiles
 pos.plus.1 <- sort(sample(1:100,number.plus.1,replace = F))
 
 acum <- c(0)
 for(i in 1:number.tiles)
 {
  if(i %in% pos.plus.1)
  {
   acum <- c(acum,acum[i]+step.len+1)
  } else
  {
   acum <- c(acum,acum[i]+step.len)
  }
 }
 
 acum[1] <- 1
 
 n.tile.signal.vector <- vector(mode = "numeric",length = number.tiles)
 for(i in 2:length(acum))
 {
  n.tile.signal.vector[(i-1)] <- mean(signal.vector[acum[i-1]:acum[i]])
 }
 
 return(n.tile.signal.vector)
}

metageneplot <- function(gene.names,genes,bw.file.name,ntile.txs,ntile.gene.body)
{
 res <- matrix(nrow=length(gene.names),ncol=2*ntile.txs+ntile.gene.body)
 rownames(res) <- gene.names
 
 bw.file <- import.bw(con = bw.file.name)
 bw.file <- as.data.frame(bw.file)
 
 for(i in 1:length(gene.names))
 {
  print(i)
  current.gene.chr <- genes[gene.names[i],"seqnames"]
  current.gene.start <- genes[gene.names[i],"start"]
  current.gene.end <- genes[gene.names[i],"end"]
  current.gene.strand <- genes[gene.names[i],"strand"]
  
  gene.body.signal <- subset(bw.file, seqnames == current.gene.chr & 
                              start >= current.gene.start & 
                              end <= current.gene.end)[["score"]]

  upstream.start.signal <- subset(bw.file, seqnames == current.gene.chr & 
                                   start >= (current.gene.start - 2000) & 
                                   end <= current.gene.start)[["score"]]
  
  downstream.end.signal <- subset(bw.file, seqnames == current.gene.chr & 
                                   start >= current.gene.end & 
                                   end <= (current.gene.end + 2000))[["score"]]
  
  ntile.signal.complete <- 
   c(ntile.signal(signal.vector = upstream.start.signal,number.tiles = ntile.txs),
     ntile.signal(signal.vector = gene.body.signal,number.tiles = ntile.gene.body),
     ntile.signal(signal.vector = downstream.end.signal,number.tiles = ntile.txs))
  
  if(current.gene.strand == "+")
  {
   res[i,] <- ntile.signal.complete
  } else
  {
   res[i,] <- rev(ntile.signal.complete)
  }
 }
 
 return(res)
}

i <- 1

gene.names[i]
plot(res[i,],type="l")
i <- i + 1


res1 <- metageneplot(gene.names = k27.gene.body[-(44:75)],genes = ostta.genes,bw.file.name = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",ntile.txs = 20,ntile.gene.body = 50)
res1 <- metageneplot(gene.names = k27.gene.body,genes = ostta.genes,bw.file.name = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",ntile.txs = 10,ntile.gene.body = 20)
plot(colMeans(res1,na.rm = T),type="l",lwd=3,col="blue",ylim=c(0,60))

plot(apply(X = res1,MARGIN = 2,FUN = median,na.rm = T),type="l",lwd=3,col="blue",ylim=c(0,60))



lines(colMeans(res1,na.rm = T),type="l",lwd=3,col="blue")#,ylim=c(0,60))




k27.internal.gene.body <- (unlist(sapply(X = subset(ll.20.peak.annotation, peak.annotation == "internal.gene.body")[["target.genes"]],
                                FUN = function(x) {strsplit(x,",")[[1]]})))
names(k27.internal.gene.body) <- NULL
length(k27.internal.gene.body)

res2 <- metageneplot(gene.names = k27.internal.gene.body,genes = ostta.genes,bw.file.name = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",ntile.txs = 20,ntile.gene.body = 50)
lines(colMeans(res2,na.rm = T),type="l",lwd=3,col="red")
lines(apply(X = res2,MARGIN = 2,FUN = median,na.rm = T),type="l",lwd=3,col="red")



k27.tss <- (unlist(sapply(X = subset(ll.20.peak.annotation, peak.annotation == "TSS")[["target.genes"]],
                                         FUN = function(x) {strsplit(x,",")[[1]]})))
names(k27.tss) <- NULL
length(k27.tss)

res3 <- metageneplot(gene.names = k27.tss,genes = ostta.genes,bw.file.name = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",ntile.txs = 20,ntile.gene.body = 50)
lines(colMeans(res3,na.rm = T),type="l",lwd=3,col="green")
lines(apply(X = res3,MARGIN = 2,FUN = median,na.rm = T),type="l",lwd=3,col="green")


k27.tes <- (unlist(sapply(X = subset(ll.20.peak.annotation, peak.annotation == "TES")[["target.genes"]],
                          FUN = function(x) {strsplit(x,",")[[1]]})))
names(k27.tes) <- NULL
length(k27.tes)

res4 <- metageneplot(gene.names = k27.tes,genes = ostta.genes,bw.file.name = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",ntile.txs = 20,ntile.gene.body = 50)
lines(colMeans(res4,na.rm = T),type="l",lwd=3,col="black")
lines(apply(X = res4,MARGIN = 2,FUN = median,na.rm = T),type="l",lwd=3,col="black")



smooth.res1 <- smooth.spline(x=1:40,y=apply(X = res1,MARGIN = 2,FUN = median,na.rm = T))
smooth.res2 <- smooth.spline(x=1:90,y=apply(X = res2,MARGIN = 2,FUN = median,na.rm = T))
smooth.res3 <- smooth.spline(x=1:90,y=apply(X = res3,MARGIN = 2,FUN = median,na.rm = T))
smooth.res4 <- smooth.spline(x=1:90,y=apply(X = res4,MARGIN = 2,FUN = median,na.rm = T))


smooth.res1$x
plot(smooth.res1$y,type="l")


res[,12]
ntile.signal.complete <- 
 c(ntile.signal(signal.vector = upstream.start.signal,number.tiles = 10),
   ntile.signal(signal.vector = gene.body.signal,number.tiles = 50),
   ntile.signal(signal.vector = downstream.end.signal,number.tiles = 10))
length(ntile.signal.complete)

plot(ntile.signal.complete,type="l")

if(current.gene.strand == "+")
{
 current.gene.tss <- current.gene.start
 current.gene.tes <- current.gene.end
} else
{
 rev(ntile.signal.complete)
}


intergenic.lens <- vector(mode = "numeric",length = nrow(ostta.genes)-1)
for(i in 1:(nrow(ostta.genes)-1))
{
 intergenic.lens[i] <- ostta.genes[(i+1),"start"] - ostta.genes[i,"end"]
}

intergenic.lens[intergenic.lens < 0] <- - intergenic.lens[intergenic.lens < 0] 
sum(intergenic.lens < 0)
mean(intergenic.lens)
summary(intergenic.lens)

plot(smooth.spline(x=1:40,y=median.metageneplot.data.k27.gene.body,spar = 0.5),
     type="l",col="blue",lwd=3,ylim=c(0,45),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)

plot(smooth.spline(x=1:40,y=median.metageneplot.data.k27.internal.gene.body),
     type="l",col="cyan",lwd=3,ylim=c(0,10),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)

plot(median.metageneplot.data.k27.internal.gene.body.3,
     type="l",col="cyan",lwd=3,ylim=c(0,10),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)



plot(median.metageneplot.data.k27.tss.3,
     type="l",col="green",lwd=3,ylim=c(0,20),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)

plot(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tss,spar = 0.5),
     type="l",col="green",lwd=3,ylim=c(0,20),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)

res <- smooth.spline(x=1:40,y=median.metageneplot.data.k27.tss)

res$y

jpeg(filename = "images/h3k27me3_metageneplots.jpg",res = 250,width = 900,height = 900)
plot(smooth.spline(x=1:40,y=median.metageneplot.data.k27.gene.body,spar = 0.5),
     type="l",col="blue",lwd=3,ylim=c(0,45),axes=F,xlab="",ylab="")
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.internal.gene.body,spar = 0.5),
      type="l",col="cyan",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tss,spar = 0.5),
      type="l",col="green",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tes,spar = 0.5),
      type="l",col="red",lwd=3)
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)
dev.off()


max.metageneplot.data.k27.gene.body.1 <- 
 apply(X = metageneplot.data.k27.gene.body.1,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.gene.body.2 <- 
 apply(X = metageneplot.data.k27.gene.body.2,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.gene.body.3 <- 
 apply(X = metageneplot.data.k27.gene.body.3,
       MARGIN = 1,
       FUN = max,na.rm = T)

max.metageneplot.data.k27.gene.body <- (max.metageneplot.data.k27.gene.body.1 + 
                                   max.metageneplot.data.k27.gene.body.2 +
                                   max.metageneplot.data.k27.gene.body.3)/3

max.metageneplot.data.k27.internal.gene.body.1 <- 
 apply(X = metageneplot.data.k27.internal.gene.body.1,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.internal.gene.body.2 <- 
 apply(X = metageneplot.data.k27.internal.gene.body.2,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.internal.gene.body.3 <- 
 apply(X = metageneplot.data.k27.internal.gene.body.3,
       MARGIN = 1,
       FUN = max,na.rm = T)

max.metageneplot.data.k27.internal.gene.body <- (max.metageneplot.data.k27.internal.gene.body.1 + 
                                         max.metageneplot.data.k27.internal.gene.body.2 +
                                         max.metageneplot.data.k27.internal.gene.body.3)/3




max.metageneplot.data.k27.tss.1 <- 
 apply(X = metageneplot.data.k27.tss.1,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.tss.2 <- 
 apply(X = metageneplot.data.k27.tss.2,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.tss.3 <- 
 apply(X = metageneplot.data.k27.tss.3,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.tss <- (max.metageneplot.data.k27.tss.1 + 
                                   max.metageneplot.data.k27.tss.2 +
                                   max.metageneplot.data.k27.tss.3)/3

max.metageneplot.data.k27.tes.1 <- 
 apply(X = metageneplot.data.k27.tes.1,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.tes.2 <- 
 apply(X = metageneplot.data.k27.tes.2,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.tes.3 <- 
 apply(X = metageneplot.data.k27.tes.3,
       MARGIN = 1,
       FUN = max,na.rm = T)
max.metageneplot.data.k27.tes <- (max.metageneplot.data.k27.tes.1 + 
                                  max.metageneplot.data.k27.tes.2 +
                                  max.metageneplot.data.k27.tes.3)/3


shapiro.test(max.metageneplot.data.k27.gene.body)
shapiro.test(max.metageneplot.data.k27.internal.gene.body)
shapiro.test(max.metageneplot.data.k27.tss)
shapiro.test(max.metageneplot.data.k27.tes)

wilcox.test(x = max.metageneplot.data.k27.gene.body,y = max.metageneplot.data.k27.internal.gene.body,alternative = "greater")$p.value
wilcox.test(x = max.metageneplot.data.k27.gene.body,y = max.metageneplot.data.k27.tss,alternative = "greater")
wilcox.test(x = max.metageneplot.data.k27.gene.body,y = max.metageneplot.data.k27.tes,alternative = "greater")
wilcox.test(x = max.metageneplot.data.k27.tss,y = max.metageneplot.data.k27.tes)
wilcox.test(x = max.metageneplot.data.k27.tss,y = max.metageneplot.data.k27.tes)


non.k27.genes <- setdiff(rownames(gene.expression),c(k27.gene.body,k27.tss,k27.tes,k27.internal.gene.body))

jpeg(filename = "images/h3k27me3_expression_barplot_different_genes.jpg",res = 250,width = 600,height = 1000)
boxplot(gene.expression.20[k27.gene.body],
        gene.expression.20[k27.tss],
        gene.expression.20[k27.internal.gene.body],
        gene.expression.20[k27.tes],
        gene.expression.20[non.k27.genes],
        col=c("blue","green","cyan","red","grey"),
        outline=F)
dev.off()

wilcox.test(x = gene.expression.20[k27.gene.body], 
            y = gene.expression.20[k27.tss],
            alternative = "less")
wilcox.test(x = gene.expression.20[k27.tss], 
            y = gene.expression.20[k27.internal.gene.body],
            alternative = "less")

wilcox.test(x = gene.expression.20[k27.internal.gene.body],
            y = gene.expression.20[k27.tes],
            alternative = "less")

wilcox.test(x = gene.expression.20[k27.tes],
            y = gene.expression.20[non.k27.genes],
            alternative = "less")


all.k27.genes

library(clusterProfiler)
library(enrichplot)
library(org.Otauriv5.eg.db)
enrichment <- enrichGO(gene = all.k27.genes, pvalueCutoff = 1,qvalueCutoff = 1,
                       OrgDb = org.Otauriv5.eg.db,
                       ont = "BP", 
                       keyType = "GID")
enrichment.df <- as.data.frame(enrichment)

write.table(x = enrichment.df, file = "go_enrichment_k27.tsv",quote = F,sep = "\t",row.names = F)

enrichment.df[which.max(enrichment.df$Count),]

sort(enrichment.df$Count,decreasing = T)[10:20]

enrichment.df[enrichment.df$Count == 26,]
enrichment.df[enrichment.df$Count == 24,]
enrichment.df[enrichment.df$Count == 22,]
enrichment.df[enrichment.df$Count == 21,]
enrichment.df[enrichment.df$Count == 20,]
enrichment.df[enrichment.df$Count == 18,]

jpeg(filename = "images/h3k27me3_treeplot.jpg",res = 250,width = 3600,height = 2000)
treeplot(x = pairwise_termsim(enrichment))
dev.off()


jpeg(filename = "images/h3k27me3_emapplot.jpg",res = 250,width = 2000,height = 2000)
emapplot(pairwise_termsim(enrichment))
dev.off()


cnetplot(pairwise_termsim(enrichment))

barplot(enrichment)


library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000003","reproduction",1.1154428097959028,0.35455656429893245,1,-0,"reproduction"),
                     c("GO:0001887","selenium compound metabolic process",0.0031524026730797234,0.33700103077893184,0.9533150770650385,0.05458229,"selenium compound metabolic process"),
                     c("GO:0002376","immune system process",0.9427113541805001,0.3647586988245658,1,-0,"immune system process"),
                     c("GO:0006099","tricarboxylic acid cycle",0.5040640114710979,0.40392617241762396,0.926776665398285,0.08410724,"tricarboxylic acid cycle"),
                     c("GO:0006103","2-oxoglutarate metabolic process",0.017201421622692253,1.0244715302516616,0.900036020673125,0.06139517,"2-oxoglutarate metabolic process"),
                     c("GO:0006525","arginine metabolic process",0.425056765430777,0.6822743332528265,0.8414193966872098,0.31651609,"2-oxoglutarate metabolic process"),
                     c("GO:0006526","arginine biosynthetic process",0.2840884163423315,0.6822743332528265,0.8197776139390401,0.55006631,"2-oxoglutarate metabolic process"),
                     c("GO:0006631","fatty acid metabolic process",1.936945636803579,0.661709121736741,0.8363267566648038,0.47549756,"2-oxoglutarate metabolic process"),
                     c("GO:0019344","cysteine biosynthetic process",0.13937119746048948,0.3647586988245658,0.8182547791420745,0.64560862,"2-oxoglutarate metabolic process"),
                     c("GO:0019541","propionate metabolic process",0.06516280052453614,0.3121779710446768,0.8755895054022355,0.56447247,"2-oxoglutarate metabolic process"),
                     c("GO:0030258","lipid modification",0.7621839254181456,0.3121779710446768,0.9072895155224819,0.63919224,"2-oxoglutarate metabolic process"),
                     c("GO:0042128","nitrate assimilation",0.06915814433459262,0.39602842215015766,0.872968299057155,0.3471953,"2-oxoglutarate metabolic process"),
                     c("GO:0043648","dicarboxylic acid metabolic process",1.2968131796952338,0.34085777107681914,0.8591201035057788,0.53948233,"2-oxoglutarate metabolic process"),
                     c("GO:0046487","glyoxylate metabolic process",0.1082932886998466,0.42431833266535085,0.8788025002701066,0.59108986,"2-oxoglutarate metabolic process"),
                     c("GO:0007017","microtubule-based process",1.4522840599239462,2.9055417594493695,0.9846414615040694,0.01131391,"microtubule-based process"),
                     c("GO:0007018","microtubule-based movement",0.6202173565622278,3.3618830741865118,0.9760729714791799,0,"microtubule-based movement"),
                     c("GO:0007049","cell cycle",2.8073119376140743,1.055588520937294,0.9834821578378161,0.01350265,"cell cycle"),
                     c("GO:0007154","cell communication",9.671968224175414,0.47456381848456913,0.9807530444684127,0.0178355,"cell communication"),
                     c("GO:0009267","cellular response to starvation",0.27913921808780134,0.9270279660356058,0.8817472625874171,0.00961211,"cellular response to starvation"),
                     c("GO:0000160","phosphorelay signal transduction system",2.0739975137672992,0.6100474021856732,0.8687040167686203,0.43251558,"cellular response to starvation"),
                     c("GO:0007166","cell surface receptor signaling pathway",1.7699816731780171,0.5778037246528189,0.8709358661055234,0.56185435,"cellular response to starvation"),
                     c("GO:0071496","cellular response to external stimulus",0.6305470826004429,0.7229162230356428,0.9204539188066697,0.65860564,"cellular response to starvation"),
                     c("GO:0015074","DNA integration",0.7252350055806647,2.657762206770456,0.8035461550887487,0.01150715,"DNA integration"),
                     c("GO:0000271","polysaccharide biosynthetic process",0.8445801891875397,0.42431833266535085,0.8524365468053611,0.2912699,"DNA integration"),
                     c("GO:0005975","carbohydrate metabolic process",5.3400641443698875,0.6489525204320686,0.9194445968125812,0.10937561,"DNA integration"),
                     c("GO:0006259","DNA metabolic process",5.6268046211314795,2.0741453286090206,0.7883371559851247,0.44302427,"DNA integration"),
                     c("GO:0006270","DNA replication initiation",0.18677061560434888,1.524614033979224,0.8223941766176285,0.48453834,"DNA integration"),
                     c("GO:0006281","DNA repair",2.895534847137308,0.7716708877524712,0.7159639044019359,0.64562006,"DNA integration"),
                     c("GO:0006298","mismatch repair",0.21782141222315915,0.6187036479778791,0.7758261295536709,0.68512103,"DNA integration"),
                     c("GO:0006468","protein phosphorylation",1.0150465485898383,2.0478177071574097,0.8059994933400036,0.18117637,"DNA integration"),
                     c("GO:0006486","protein glycosylation",0.7526798873357411,0.6623793059324244,0.7700801528471054,0.66477721,"DNA integration"),
                     c("GO:0006508","proteolysis",5.2622572267907,2.215841795989129,0.8331839915576025,0.43614593,"DNA integration"),
                     c("GO:0006637","acyl-CoA metabolic process",0.5539971825062208,0.39602842215015766,0.766633541793333,0.51243251,"DNA integration"),
                     c("GO:0008213","protein alkylation",0.000579213939150692,0.33320452706495884,0.8948472249778151,0.37850183,"DNA integration"),
                     c("GO:0009100","glycoprotein metabolic process",0.9570906483124367,0.6014649273302938,0.8398249210805684,0.35568771,"DNA integration"),
                     c("GO:0009119","ribonucleoside metabolic process",0.3506733719429951,0.39602842215015766,0.7958005207549018,0.45483879,"DNA integration"),
                     c("GO:0009226","nucleotide-sugar biosynthetic process",0.35550672740969513,0.3033112734255381,0.7919369046997367,0.48997503,"DNA integration"),
                     c("GO:0010501","RNA secondary structure unwinding",0.0022996025754365773,0.33700103077893184,0.8841308449932884,0.22153339,"DNA integration"),
                     c("GO:0016266","O-glycan processing",0.009247705105078281,0.6822743332528265,0.832553230099121,0.68426048,"DNA integration"),
                     c("GO:0016310","phosphorylation",5.235381700014107,0.9200287559031499,0.8995173712742447,0.51243059,"DNA integration"),
                     c("GO:0016579","protein deubiquitination",0.0016070105886223452,1.363169883748535,0.8876399511853289,0.40321557,"DNA integration"),
                     c("GO:0018202","peptidyl-histidine modification",1.0281836136778093,0.9270279660356058,0.8252111825191477,0.68739488,"DNA integration"),
                     c("GO:0018342","protein prenylation",0.00018239077232830298,0.3647586988245658,0.8956376571084298,0.35393725,"DNA integration"),
                     c("GO:0019439","aromatic compound catabolic process",2.2395984760807317,0.488426217953643,0.8298501364385389,0.67803842,"DNA integration"),
                     c("GO:0032543","mitochondrial translation",0.2221051306327077,0.7671257790005732,0.8301308068014243,0.30712983,"DNA integration"),
                     c("GO:0034404","nucleobase-containing small molecule biosynthetic process",0.2580336480412168,0.3647586988245658,0.8075679421637405,0.47645433,"DNA integration"),
                     c("GO:0034655","nucleobase-containing compound catabolic process",1.6564680648053147,0.5761569684631386,0.8082436737396426,0.50011278,"DNA integration"),
                     c("GO:0035383","thioester metabolic process",0.5539971825062208,0.39602842215015766,0.9156152717960674,0.69278956,"DNA integration"),
                     c("GO:0042398","cellular modified amino acid biosynthetic process",0.5343877097407617,0.40392617241762396,0.8708997557571673,0.19203267,"DNA integration"),
                     c("GO:0044270","cellular nitrogen compound catabolic process",2.004743244566209,0.4588180186234574,0.8286580702160516,0.69482359,"DNA integration"),
                     c("GO:0046434","organophosphate catabolic process",0.9710632603168847,0.868457312023167,0.8769887779688311,0.41091092,"DNA integration"),
                     c("GO:0046700","heterocycle catabolic process",2.042108705491845,0.4685002750565224,0.8310802481954842,0.69649215,"DNA integration"),
                     c("GO:0070646","protein modification by small protein removal",0.0641152859537317,1.1763933331577958,0.8585287404338742,0.52767687,"DNA integration"),
                     c("GO:0071897","DNA biosynthetic process",0.7667855953091854,2.0080982132635268,0.7764795050343053,0.55438786,"DNA integration"),
                     c("GO:0090304","nucleic acid metabolic process",13.541109943481535,0.5761569684631386,0.7844903754380635,0.55585571,"DNA integration"),
                     c("GO:0140053","mitochondrial gene expression",0.3272164397774754,0.7229162230356428,0.8689292553283973,0.16191549,"DNA integration"),
                     c("GO:1901361","organic cyclic compound catabolic process",2.407836710372425,0.4588180186234574,0.8509596622336734,0.55813667,"DNA integration"),
                     c("GO:0016441","regulation of gene expression",0.049918382594123255,0.7229162230356428,0.9152721684520655,0.09459438,"regulation of gene expression"),
                     c("GO:0006450","regulation of translational fidelity",0.29561600610151356,0.3647586988245658,0.9565840743545636,0.13549893,"regulation of gene expression"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6811968028967992,0.5129770046876073,0.9187967357951126,0.29038639,"regulation of gene expression"),
                     c("GO:0010629","negative regulation of gene expression",1.1219250764340574,0.531092579249663,0.891938194299043,0.58730682,"regulation of gene expression"),
                     c("GO:0017148","negative regulation of translation",0.21779183534115995,0.33700103077893184,0.9025830172260265,0.66877869,"regulation of gene expression"),
                     c("GO:0031326","regulation of cellular biosynthetic process",14.068660464000152,0.3894214410778525,0.8933198148321295,0.55551169,"regulation of gene expression"),
                     c("GO:0033044","regulation of chromosome organization",0.18193972487781543,0.3121779710446768,0.9545386992795665,0.15097414,"regulation of gene expression"),
                     c("GO:0034249","negative regulation of amide metabolic process",0.23996217313971518,0.33700103077893184,0.9072896668585831,0.68743154,"regulation of gene expression"),
                     c("GO:1903311","regulation of mRNA metabolic process",0.4779772015478174,0.3647586988245658,0.9285627652482316,0.34495101,"regulation of gene expression"),
                     c("GO:2000241","regulation of reproductive process",0.1639520511419782,0.3121779710446768,0.9569481717332748,0.14866641,"regulation of gene expression"),
                     c("GO:0030029","actin filament-based process",0.8246601591610891,0.6822743332528265,0.9855133147715953,0.01167087,"actin filament-based process"),
                     c("GO:0030154","cell differentiation",2.279586420543629,0.4717135495803034,0.9548569791025598,0.01433045,"cell differentiation"),
                     c("GO:0009888","tissue development",0.6842315881691289,0.3494179598129715,0.9723962948998165,0.69552749,"cell differentiation"),
                     c("GO:0048229","gametophyte development",0.017253181166190824,0.3647586988245658,0.9784815653260818,0.48986076,"cell differentiation"),
                     c("GO:0048507","meristem development",0.03934957675974807,0.33700103077893184,0.9773200335076074,0.52461364,"cell differentiation"),
                     c("GO:0040011","locomotion",0.5148398554794674,2.9956162519741056,1,-0,"locomotion"),
                     c("GO:0043547","positive regulation of GTPase activity",0.0008897712001421268,0.9270279660356058,0.945105760674195,-0,"positive regulation of GTPase activity"),
                     c("GO:0051168","nuclear export",0.36225518598584233,0.8822380235820269,0.9193145237527721,0.00984617,"nuclear export"),
                     c("GO:0006403","RNA localization",0.38503431460554743,0.5719097997394266,0.9492853307589769,0.48954515,"nuclear export"),
                     c("GO:0006614","SRP-dependent cotranslational protein targeting to membrane",0.1731775735855571,0.3647586988245658,0.9239488841068346,0.45707974,"nuclear export"),
                     c("GO:0006890","retrograde vesicle-mediated transport, Golgi to endoplasmic reticulum",0.1396694310206479,0.3121779710446768,0.9533956800814197,0.22525423,"nuclear export"),
                     c("GO:0015931","nucleobase-containing compound transport",0.5729682875685309,0.3245363289714383,0.953096561853417,0.48659654,"nuclear export"),
                     c("GO:0031503","protein-containing complex localization",0.2928998624379218,0.3647586988245658,0.9641396486534709,0.22721235,"nuclear export"),
                     c("GO:0033750","ribosome localization",0.08453565823400183,0.39602842215015766,0.9615600081516522,0.20497386,"nuclear export"),
                     c("GO:0050658","RNA transport",0.3380539022900098,0.6187036479778791,0.9329722181898726,0.24285839,"nuclear export"),
                     c("GO:0051169","nuclear transport",0.5733749696960196,0.5066832555406448,0.9263142426302791,0.63305059,"nuclear export"),
                     c("GO:0051301","cell division",1.5693197819947182,0.5659276631444706,0.9845137269828846,0.01363807,"cell division"),
                     c("GO:0061982","meiosis I cell cycle process",0.22371953544182982,1.4436467858208988,0.8829935432772962,0.00942192,"meiosis I cell cycle process"),
                     c("GO:0098813","nuclear chromosome segregation",0.33031461816688995,1.0576395600847068,0.9027750717677788,0.64747101,"meiosis I cell cycle process"),
                     c("GO:0070085","glycosylation",0.8372204750500761,0.5642517085182833,0.9603127258883271,0.04818597,"glycosylation"),
                     c("GO:0071103","plasma membrane bounded cell projection organization",0.8233982121957907,2.48892688143064,0.8072006030847753,0.01166889,"plasma membrane bounded cell projection organization"),
                     c("GO:0000027","ribosomal large subunit assembly",0.012821578346646382,0.8154250425427645,0.8458879193561408,0.48240062,"plasma membrane bounded cell projection organization"),
                     c("GO:0000245","spliceosomal complex assembly",0.1148766096848317,0.33700103077893184,0.7176266579933399,0.64218601,"plasma membrane bounded cell projection organization"),
                     c("GO:0000280","nuclear division",0.45164391760787703,0.4717135495803034,0.8236095327076562,0.491872,"plasma membrane bounded cell projection organization"),
                     c("GO:0006415","translational termination",0.15480047090339727,0.3647586988245658,0.7363259627235162,0.48918433,"plasma membrane bounded cell projection organization"),
                     c("GO:0006479","protein methylation",0.000579213939150692,0.33320452706495884,0.8923844376711241,0.51387372,"plasma membrane bounded cell projection organization"),
                     c("GO:0006996","organelle organization",6.781836087483291,0.901501456180904,0.8241728840064373,0.52336479,"plasma membrane bounded cell projection organization"),
                     c("GO:0007010","cytoskeleton organization",1.98151060375585,1.5080424402092374,0.8150504345607171,0.61927091,"plasma membrane bounded cell projection organization"),
                     c("GO:0030030","cell projection organization",1.3282977705833654,1.6133625971416043,0.8512620269157171,0.41569968,"plasma membrane bounded cell projection organization"),
                     c("GO:0042255","ribosome assembly",0.2624923630025938,0.42431833266535085,0.8171364501582283,0.65665084,"plasma membrane bounded cell projection organization"),
                     c("GO:0043248","proteasome assembly",0.04973352708162835,0.39602842215015766,0.8611074492352145,0.52931974,"plasma membrane bounded cell projection organization"),
                     c("GO:0048285","organelle fission",0.5487596096521986,0.30148135432836,0.8356875754692098,0.50144746,"plasma membrane bounded cell projection organization"),
                     c("GO:0051258","protein polymerization",0.1780454354145895,1.0725168728382306,0.8473253852004743,0.32972263,"plasma membrane bounded cell projection organization"),
                     c("GO:0051276","chromosome organization",1.5299726699751368,1.9364166059501378,0.8195931631949569,0.55870698,"plasma membrane bounded cell projection organization"),
                     c("GO:0070475","rRNA base methylation",0.22436283262531206,0.33700103077893184,0.6953574452494768,0.4998146,"plasma membrane bounded cell projection organization"),
                     c("GO:0070897","transcription preinitiation complex assembly",0.11026754557329214,0.39602842215015766,0.7195583694392909,0.56139216,"plasma membrane bounded cell projection organization"),
                     c("GO:0070925","organelle assembly",1.3561320812847646,0.6413546437511954,0.8010719571331075,0.5915788,"plasma membrane bounded cell projection organization"),
                     c("GO:0097435","supramolecular fiber organization",0.801568008540226,0.5477838551470221,0.8581207014692784,0.38125869,"plasma membrane bounded cell projection organization"),
                     c("GO:0120036","plasma membrane bounded cell projection organization",1.1030821378604103,1.7112755303804414,0.8324651987784939,0.3943391,"plasma membrane bounded cell projection organization"),
                     c("GO:0097354","prenylation",0.009942761832059112,0.3647586988245658,0.9497005745047677,0.05901557,"prenylation"),
                     c("GO:2001057","reactive nitrogen species metabolic process",0.11012952012396263,0.33700103077893184,0.9373241093491478,0.09288199,"reactive nitrogen species metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
#pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches
jpeg(filename = "images/h3k27me3_revigo_treemap.jpg",res = 300,width = 1400,height = 900)

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
 stuff,
 index = c("representative","description"),
 vSize = "value",
 type = "categorical",
 vColor = "representative",
 title = "",
 inflate.labels = TRUE,      # set this to TRUE for space-filling group labels - good for posters
 lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
 #bg.labels = "#CCCCCCAA",   # define background color of group labels
 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
 position.legend = "none"
)

dev.off()



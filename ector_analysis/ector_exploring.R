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


library(rtracklayer)
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

enrichment.df[enrichment.df$Description == "cilium organization",]
enrichment.df[enrichment.df$Description == "locomotion",]


enrichment.df[enrichment.df$Description == "regulation of gene expression",]


barplot(c(gene.expression.20["ostta14g00070"],gene.expression.10["ostta14g00070"]))
barplot(c(gene.expression.20["ostta06g04460"],gene.expression.10["ostta06g04460"]))

gene.id <- "ostta11g02490"
barplot(c(gene.expression.20[gene.id],gene.expression.10[gene.id]))

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



k4.ll.20
ll.20

enrichment <- enrichGO(gene = unique(c(k27.gene.body,k27.tss)), pvalueCutoff = 1,qvalueCutoff = 1,
                       OrgDb = org.Otauriv5.eg.db,
                       ont = "BP", 
                       keyType = "GID")
enrichment.df <- as.data.frame(enrichment)


treeplot(x = pairwise_termsim(enrichment))


gene.expression.10 <- 
 (temp.gene.expression$T10_1 + 
   temp.gene.expression$T10_2 + 
   temp.gene.expression$T10_3)/3

gene.expression.20[""]

count.1 <- read.table(file = "chipseq_data/bed_data/h3k27me3/temp20c/k27_20_ll_chip_1_counts.bed",header = F,sep = "\t")
count.2 <- read.table(file="chipseq_data/bed_data/h3k27me3/temp20c/k27_20_ll_chip_2n_counts.bed",header = F,sep = "\t")
count.3 <- read.table(file="chipseq_data/bed_data/h3k27me3/temp20c/k27_20_ll_chip_3n_counts.bed",header = F,sep = "\t")
count.1 <- count.1[,c(1:3,10)]
count.2 <- count.2[,c(1:3,10)]
count.3 <- count.3[,c(1:3,10)]

mapped.reads.k27.20c.1 <- 3146662/1E6 
mapped.reads.k27.20c.2 <- 5697849/1E6
mapped.reads.k27.20c.3 <- 5157582/1E6

rpm.1 <- count.1
rpm.2 <- count.2
rpm.3 <- count.3

rpm.1[,4] <- rpm.1[,4]/mapped.reads.k27.20c.1
rpm.2[,4] <- rpm.2[,4]/mapped.reads.k27.20c.2
rpm.3[,4] <- rpm.3[,4]/mapped.reads.k27.20c.3


plot(log2(rpm.1[,4]),log2(rpm.2[,4]))
plot(log2(rpm.1[,4]),log2(rpm.3[,4]))
plot(log2(rpm.2[,4]),log2(rpm.3[,4]))

rpm <- rpm.1
rpm[,4] <- (rpm.1[,4] + rpm.2[,4] + rpm.3[,4])/3

rpkm <- rpm
rpkm[,4] <- (rpm[,4]/(rpm[,3]-rpm[,2]))

summary(rpkm[,4])
quantile(rpkm[,4])
genes.020 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.20)]
genes.020 <- genes.020[genes.020 != ""]

genes.030 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.30)]
genes.030 <- genes.030[genes.030 != ""]
genes.030 <- unlist(strsplit(genes.030,split = ","))
length(genes.030)

genes.060 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.60 & rpkm[,4] > 0.30)]
genes.060 <- genes.060[genes.060 != ""]
genes.060 <- unlist(strsplit(genes.060,split = ","))
length(genes.060)

genes.090 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] > 0.60)]
genes.090 <- genes.090[genes.090 != ""]
genes.090 <- unlist(strsplit(genes.090,split = ","))
length(genes.090)


boxplot(gene.expression.20[genes.030],
        gene.expression.20[genes.060],
        gene.expression.20[genes.090],outline=F)

jpeg(filename = "../images/boxplot_chip_signal_gene_expression.jpg",width = 400,height = 600)
boxplot(gene.expression.20[genes.2nd.quant],
        gene.expression.20[genes.3rd.quant],
        gene.expression.20[genes.4th.quant],
        names=c("RPKM < 0.3","0.3 < RPKM < 0.5", "RPKM > 0.5"),
        ylab="TPM", cex.lab=1.5, outline=F, col=c("lightblue","blue","darkblue"))
dev.off()


genes.020 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.20)]
genes.020 <- genes.020[genes.020 != ""]
genes.020 <- unlist(strsplit(genes.020,split = ","))
length(genes.020)

genes.030 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.30)]
genes.030 <- genes.030[genes.030 != ""]
genes.030 <- unlist(strsplit(genes.030,split = ","))
length(genes.030)

genes.040 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.40 & rpkm[,4] >= 0.30)]
genes.040 <- genes.040[genes.040 != ""]
genes.040 <- unlist(strsplit(genes.040,split = ","))
length(genes.040)


genes.050 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.50 & rpkm[,4] >= 0.40)]
genes.050 <- genes.050[genes.050 != ""]
genes.050 <- unlist(strsplit(genes.050,split = ","))
length(genes.050)

genes.060 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] >= 0.60)]
genes.060 <- genes.060[genes.060 != ""]
genes.060 <- unlist(strsplit(genes.060,split = ","))
length(genes.060)


genes.m060 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] >= 0.60)]
genes.m060 <- genes.m060[genes.m060 != ""]
genes.m060 <- unlist(strsplit(genes.m060,split = ","))
length(genes.m060)



boxplot(gene.expression.20[genes.030],
        gene.expression.20[genes.040],
        gene.expression.20[genes.050],
        gene.expression.20[genes.060]
        #gene.expression.20[genes.m060]
        ,outline=F)

quantile(rpkm[ll.20.peak.annotation$peak.annotation == "gene.body",4])

genes.040 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.40)]
genes.040 <- genes.040[genes.040 != ""]
genes.040 <- unlist(strsplit(genes.040,split = ","))
genes.040 <- intersect(genes.040, k27.gene.body)
length(genes.040)

genes.050 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.60 & rpkm[,4] >= 0.40)]
genes.050 <- genes.050[genes.050 != ""]
genes.050 <- unlist(strsplit(genes.050,split = ","))
genes.050 <- intersect(genes.050, k27.gene.body)
length(genes.050)

genes.060 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.80 & rpkm[,4] >= 0.60)]
genes.060 <- genes.060[genes.060 != ""]
genes.060 <- unlist(strsplit(genes.060,split = ","))
genes.060 <- intersect(genes.060, k27.gene.body)
length(genes.060)

# genes.m050 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] >= 0.8)]
# genes.m050 <- genes.m050[genes.m050 != ""]
# genes.m050 <- unlist(strsplit(genes.m050,split = ","))
# genes.m050 <- intersect(genes.m050, k27.gene.body)
# length(genes.m050)

gene.set <- k27.gene.body
boxplot(gene.expression.20[intersect(genes.040,gene.set)],
        gene.expression.20[intersect(genes.050,gene.set)],
        gene.expression.20[intersect(genes.060,gene.set)]#,
#        gene.expression.20[intersect(genes.m050,gene.set)]
#        gene.expression.20[genes.m050]
        ,outline=F)


quantile(rpkm[ll.20.peak.annotation$peak.annotation == "TSS",4])

genes.1 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.25)]
genes.1 <- genes.1[genes.1 != ""]
genes.1 <- unlist(strsplit(genes.1,split = ","))
genes.1 <- intersect(genes.1, k27.tss)
length(genes.1)

genes.2 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.4 & rpkm[,4] >= 0.25)]
genes.2 <- genes.2[genes.2 != ""]
genes.2 <- unlist(strsplit(genes.2,split = ","))
genes.2 <- intersect(genes.2, k27.internal.gene.body)
length(genes.2)

genes.3 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] >= 0.4)]
genes.3 <- genes.3[genes.3 != ""]
genes.3 <- unlist(strsplit(genes.3,split = ","))
genes.3 <- intersect(genes.3, k27.tss)
length(genes.3)

gene.set <- k27.tss
boxplot(gene.expression.20[intersect(genes.1,gene.set)],
        gene.expression.20[intersect(genes.2,gene.set)],
        gene.expression.20[intersect(genes.3,gene.set)]#,
        #gene.expression.20[intersect(genes.m050,gene.set)]
        #gene.expression.20[genes.m060]
        ,outline=F)



quantile(rpkm[ll.20.peak.annotation$peak.annotation == "internal.gene.body",4])

genes.018 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.2)]
genes.018 <- genes.018[genes.018 != ""]
genes.018 <- unlist(strsplit(genes.018,split = ","))
genes.018 <- intersect(genes.018, k27.internal.gene.body)
length(genes.018)

genes.020 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.3 & rpkm[,4] >= 0.2)]
genes.020 <- genes.020[genes.020 != ""]
genes.020 <- unlist(strsplit(genes.020,split = ","))
genes.020 <- intersect(genes.020, k27.internal.gene.body)
length(genes.020)


genes.030 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.3 & rpkm[,4] >= 0.20)]
genes.030 <- genes.030[genes.030 != ""]
genes.030 <- unlist(strsplit(genes.030,split = ","))
genes.030 <- intersect(genes.030, k27.internal.gene.body)
length(genes.030)

genes.060 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] < 0.60 & rpkm[,4] >= 0.50)]
genes.060 <- genes.060[genes.060 != ""]
genes.060 <- unlist(strsplit(genes.060,split = ","))
genes.060 <- intersect(genes.060, k27.gene.body)
length(genes.060)

genes.m050 <- ll.20.peak.annotation$target.genes[which(rpkm[,4] >= 0.50)]
genes.m050 <- genes.m050[genes.m060 != ""]
genes.m050 <- unlist(strsplit(genes.m050,split = ","))
genes.m050 <- intersect(genes.m050, k27.gene.body)
length(genes.m050)


gene.set <- k27.internal.gene.body
boxplot(gene.expression.20[intersect(genes.018,gene.set)],
 gene.expression.20[intersect(genes.020,gene.set)],
        gene.expression.20[intersect(genes.030,gene.set)]#,
        #gene.expression.20[intersect(genes.m050,gene.set)]
        #gene.expression.20[genes.m060]
        ,outline=F)





quantile(rpkm[,4])


plot(log2(rpkm[,4]),log2(gene.expression.20[ll.20.peak.annotation$target.genes]+1))
plot(rpkm[,4],log2(gene.expression.20[ll.20.peak.annotation$target.genes]+1),pch=19,cex=0.5,xlim=c(0,1))
summary(lm(log2(gene.expression.20[ll.20.peak.annotation$target.genes]+1) ~ rpkm[,4]))
abline(lm(log2(gene.expression.20[ll.20.peak.annotation$target.genes]+1) ~ rpkm[,4]))



plot(rpkm[ll.20.peak.annotation$peak.annotation == "internal.gene.body",4],
     log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "internal.gene.body"]
]+1),pch=19,cex=0.5,xlim=c(0,1),main="internal")

plot(rpkm[ll.20.peak.annotation$peak.annotation == "gene.body",4],
     log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "gene.body"]
     ]+1),pch=19,cex=0.5,xlim=c(0,1), main="gene_body")


plot(rpkm[ll.20.peak.annotation$peak.annotation == "TSS",4],
     log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "TSS"]
     ]+1),pch=19,cex=0.5, main = "TSS")


library(rtracklayer)
bw.file.name.1 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw"
bw.data.1 <- import.bw(con = bw.file.name.1)
bw.data.1 <- as.data.frame(bw.data.1)

bw.file.name.2 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_2n.bw"
bw.data.2 <- import.bw(con = bw.file.name.2)
bw.data.2 <- as.data.frame(bw.data.2)

bw.file.name.3 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_3n.bw"
bw.data.3 <- import.bw(con = bw.file.name.3)
bw.data.3 <- as.data.frame(bw.data.3)


peak.signal <- vector(mode = "numeric",length = nrow(ll.20))
for(i in 1:nrow(ll.20))
{
 print(i)
 current.peak <- ll.20[i,]
 
 current.peak.chr <- current.peak[[1]]
 current.peak.start <- current.peak[[2]]
 current.peak.end <- current.peak[[3]]
 
 peak.signal.1 <- max(subset(bw.data.1, seqnames == current.peak.chr & 
                              start >= current.peak.start & 
                              end <= current.peak.end)[["score"]])
 
 peak.signal.2 <- max(subset(bw.data.2, seqnames == current.peak.chr & 
                              start >= current.peak.start & 
                              end <= current.peak.end)[["score"]])
 
 peak.signal.3 <- max(subset(bw.data.3, seqnames == current.peak.chr & 
                              start >= current.peak.start & 
                              end <= current.peak.end)[["score"]])
 
 peak.signal[i] <- mean(c(peak.signal.1,peak.signal.2,peak.signal.3))
}



plot(log2(peak.signal),log2(gene.expression.20[ll.20.peak.annotation$target.genes]+1),pch=19,cex=0.6)
summary(lm(log2(gene.expression.20[ll.20.peak.annotation$target.genes]+1) ~ peak.signal))
abline(lm(log2(gene.expression.20[ll.20.peak.annotation$target.genes]+1) ~ peak.signal))


plot(log2(peak.signal[ll.20.peak.annotation$peak.annotation == "gene.body"]),
     log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "gene.body"]
     ]+1),pch=19,cex=0.5, main="gene_body")
summary(lm(log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "gene.body"]
]+1) ~ log2(peak.signal[ll.20.peak.annotation$peak.annotation == "gene.body"])))
abline(lm(log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "gene.body"]
]+1) ~ log2(peak.signal[ll.20.peak.annotation$peak.annotation == "gene.body"])))


plot(log2(peak.signal[ll.20.peak.annotation$peak.annotation == "TSS"]),
     log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "TSS"]
     ]+1),pch=19,cex=0.5, main="TSS")
summary(lm(log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "TSS"]
]+1) ~ log2(peak.signal[ll.20.peak.annotation$peak.annotation == "TSS"])))
abline(lm(log2(gene.expression.20[ll.20.peak.annotation$target.genes[ll.20.peak.annotation$peak.annotation == "TSS"]
]+1) ~ log2(peak.signal[ll.20.peak.annotation$peak.annotation == "TSS"])))




k4.ll.20 <- read.table(file="chipseq_data/bed_data/h3k4me3/h3k4me3_ll_20.bed",
                       header = F)

k4.ll.20 <- k4.ll.20[,1:3]

nrow(k4.ll.20)

jpeg(filename = "images/h3k4me3_genome_wide_distribution_ll_20.jpg",res = 300,width = 900,height = 900)
par(mar=c(0,0,0,0))
plot(x=c(0,max.chr.len), y=c(0,4*number.chrs), col = "white", xlab = "", ylab = "", axes=F)  

for(i in 1:number.chrs)
{
 polygon(x = c(0, chr.lens[i], chr.lens[i], 0),   
         y = c(4*i+1, 4*i+1, 4*i-1, 4*i-1),    
         col = "white") 
}

for(i in 1:nrow(k4.ll.20))
{
 current.chr <- k4.ll.20[i,1]
 current.start <- k4.ll.20[i,2]
 current.end <- k4.ll.20[i,3]
 
 current.line <- (20:1)[current.chr]
 
 polygon(x = c(current.start, current.end, current.end, current.start),   
         y = c(4*current.line+1, 4*current.line+1, 4*current.line-1, 4*current.line-1),    
         col = "darkorange",border="darkorange")
}
dev.off()











h3k4me3.peaks.len <- k4.ll.20$V3 - k4.ll.20$V2
sum(h3k4me3.peaks.len)/sum(chr.lens)
k4.ll.20 <- cbind(k4.ll.20,h3k4me3.peaks.len)
colnames(k4.ll.20) <- c("chr","start","end","length")
head(k4.ll.20)
number.of.h3k4.peaks.chr <- vector(mode = "numeric",length = 20)
k4.peaks.len.chr <- vector(mode = "numeric",length = 20)

for(i in 1:20)
{
 number.of.h3k4.peaks.chr[i] <- nrow(subset(k4.ll.20, chr == i))
 k4.peaks.len.chr[i] <- sum(subset(k4.ll.20, chr == i)[,4]) 
}

percentage.chr.k4 <- k4.peaks.len.chr/rev(chr.lens)

k4.z.scores <- (percentage.chr.k4 - mean(percentage.chr.k4))/sd(percentage.chr.k4)

which(k4.z.scores > 2.5)
which(k4.z.scores < -2.5)

p.vals <- pnorm((percentage.chr.k4 - mean(percentage.chr.k4[-20]))/sd(percentage.chr.k4[-20]),
                lower.tail = F)

p.vals

values.for.colors <- round(-log10(p.vals)*40+1)
typeof(values.for.colors)
values.for.colors[20] <- 100
length(values.for.colors)
white.darkorange.palette <- colorRampPalette(c('white','darkorange'))
white.to.darkorange <- white.darkorange.palette(100)

white.to.darkorange[values.for.colors]

jpeg(filename = "images/h3k4me3_chr_barplot_ll_20.jpg",res = 300,width = 1800,height = 900)
barplot(100*percentage.chr.k4,border="black",names.arg = paste0("Chr",1:20),las=2,col=white.darkorange.palette(100)[values.for.colors]
)
abline(h = 100*mean(percentage.chr.k4[-20]),lty=5)
dev.off()




target.genes <- vector(mode = "character",length = nrow(k4.ll.20))
peak.annotation <- vector(mode = "character",length = nrow(k4.ll.20))
i <- 1
for(i in 1:nrow(k4.ll.20))
{
 current.chr <- k4.ll.20[i,1]
 current.start <- k4.ll.20[i,2]
 current.end <- k4.ll.20[i,3]
 
 complete.gene.body <- subset(ostta.genes, 
                              seqnames == current.chr & 
                               start >= current.start & 
                               end <= current.end)
 
 overlap.tss.pos.strand <- subset(ostta.genes, 
                                  seqnames == current.chr & 
                                   start >= current.start & start <= current.end & 
                                   strand == "+")
 
 overlap.tss.neg.strand <- subset(ostta.genes, 
                                  seqnames == current.chr & 
                                   end >= current.start & end <= current.end & 
                                   strand == "-")
 
 overlap.tes.pos.strand <- subset(ostta.genes, 
                                  seqnames == current.chr & 
                                   end >= current.start & end <= current.end & 
                                   strand == "+")
 
 overlap.tes.neg.strand <- subset(ostta.genes, 
                                  seqnames == current.chr & 
                                   start >= current.start & start <= current.end & 
                                   strand == "-")
 
 inside.gene.body <- subset(ostta.genes, seqnames == current.chr & 
                             start <= current.start & end >= current.end)
 
 
 if(nrow(complete.gene.body) > 0)
 {
  target.genes[i] <- paste(unlist(complete.gene.body$gene_id),collapse = ",")
  peak.annotation[i] <- "gene.body"
 } else if((nrow(overlap.tss.pos.strand) > 0) | (nrow(overlap.tss.neg.strand) > 0))
 {
  target.genes[i] <- paste(c(unlist(overlap.tss.pos.strand$gene_id),
                             unlist(overlap.tss.neg.strand$gene_id)),collapse = ",")
  peak.annotation[i] <- "TSS"
 } else if(nrow(overlap.tes.pos.strand) > 0)
 {
  target.genes[i] <- paste(unlist(overlap.tes.pos.strand$gene_id),collapse = ",")
  peak.annotation[i] <- "TES"
 } else if(nrow(overlap.tes.neg.strand) > 0)
 {
  target.genes[i] <- paste(unlist(overlap.tes.neg.strand$gene_id),collapse = ",")
  peak.annotation[i] <- "TES"
 } else if(nrow(inside.gene.body) > 0)
 {
  target.genes[i] <- paste(unlist(inside.gene.body$gene_id),collapse = ",")
  peak.annotation[i] <- "internal.gene.body"
 } else
 {
  peak.annotation[i] <- "intergenic"
 }
}

freq.peak.annotation <- table(peak.annotation)
jpeg(filename = "images/h3k4me3_pie_chart_ll_20.jpg",res = 300,width = 900,height = 900)
par(mar=c(0,0,0,0))
pie(freq.peak.annotation[c("gene.body",
                           "internal.gene.body",
                           "TSS", 
                           "TES",
                           "intergenic")],
    labels = c("Full Gene Body", 
               "Gene Body", 
               "TSS", 
               "TES", 
               "Intergenic"),
    clockwise = T,border = "black",col=c("blue","cyan","green","red","grey"))
dev.off()

percentage.freq.peak.annotation <- 100*freq.peak.annotation/sum(freq.peak.annotation)


percentage.freq.peak.annotation["gene.body"] +
 percentage.freq.peak.annotation["internal.gene.body"]

percentage.freq.peak.annotation

k4.ll.20.peak.annotation <- cbind(k4.ll.20[,1:3],data.frame(peak.annotation,target.genes))

head(k4.ll.20.peak.annotation)

write.table(x = k4.ll.20.peak.annotation,file = "tables/H3K4me3_LL_20C_peak_annotation.tsv",quote = F,sep = "\t",row.names = F)

all.k4.genes <- (unlist(sapply(X =k4.ll.20.peak.annotation$target.genes,
                                FUN = function(x) {strsplit(x,",")[[1]]})))
names(all.k4.genes) <- NULL
length(all.k4.genes)
length(all.k4.genes)/7668

k4.k27.genes <- intersect(all.k4.genes, all.k27.genes)

length(k4.k27.genes)/length(all.k4.genes)
length(k4.k27.genes)/length(all.k27.genes)






## Genes with their entire gene body marked
k4.gene.body <- (unlist(sapply(X = subset(k4.ll.20.peak.annotation, 
                                          peak.annotation == "gene.body")[["target.genes"]],
                                FUN = function(x) {strsplit(x,",")[[1]]})))
names(k4.gene.body) <- NULL
length(k4.gene.body)

## Genes with H3K27 inside their gene body with no overlap
## in the TSS or TES.
k4.internal.gene.body <- (unlist(sapply(X = subset(k4.ll.20.peak.annotation, 
                                                   peak.annotation == "internal.gene.body")[["target.genes"]],
                                         FUN = function(x) {strsplit(x,",")[[1]]})))
names(k4.internal.gene.body) <- NULL
length(k4.internal.gene.body)

## Genes with H3K27me3 overlapping their TSS
k4.tss <- (unlist(sapply(X = subset(k4.ll.20.peak.annotation, 
                                    peak.annotation == "TSS")[["target.genes"]],
                          FUN = function(x) {strsplit(x,",")[[1]]})))
names(k4.tss) <- NULL
length(k4.tss)

## Genes with H3K27me3 overlapping their TES
k4.tes <- (unlist(sapply(X = subset(k4.ll.20.peak.annotation, 
                                    peak.annotation == "TES")[["target.genes"]],
                          FUN = function(x) {strsplit(x,",")[[1]]})))
names(k4.tes) <- NULL
length(k4.tes)

metageneplot.data.k4.gene.body.1 <- 
 metageneplot(gene.names = k4.gene.body,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_1.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.gene.body.1 <- 
 apply(X = metageneplot.data.k4.gene.body.1,
       MARGIN = 2,
       FUN = median,na.rm = T)

metageneplot.data.k4.gene.body.2 <- 
 metageneplot(gene.names = k4.gene.body,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_2.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.gene.body.2 <- 
 apply(X = metageneplot.data.k4.gene.body.2,
       MARGIN = 2,
       FUN = median,na.rm = T)

metageneplot.data.k4.gene.body.3 <- 
 metageneplot(gene.names = k4.gene.body,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_3.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.gene.body.3 <- 
 apply(X = metageneplot.data.k4.gene.body.3,
       MARGIN = 2,
       FUN = median,na.rm = T)

median.metageneplot.data.k4.gene.body <- (median.metageneplot.data.k4.gene.body.1 + 
                                           median.metageneplot.data.k4.gene.body.2 + 
                                           median.metageneplot.data.k4.gene.body.3)/3

plot(smooth.spline(x=1:40,y=median.metageneplot.data.k4.gene.body,spar = 0.5),
     type="l",col="blue",lwd=3,ylim=c(0,70),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)


metageneplot.data.k4.internal.gene.body.1 <- 
 metageneplot(gene.names = k4.internal.gene.body,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_1.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.internal.gene.body.1 <- 
 apply(X = metageneplot.data.k4.internal.gene.body.1,
       MARGIN = 2,
       FUN = median,na.rm = T)


metageneplot.data.k4.internal.gene.body.2 <- 
 metageneplot(gene.names = k4.internal.gene.body,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_2.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.internal.gene.body.2 <- 
 apply(X = metageneplot.data.k4.internal.gene.body.2,
       MARGIN = 2,
       FUN = median,na.rm = T)


metageneplot.data.k4.internal.gene.body.3 <- 
 metageneplot(gene.names = k4.internal.gene.body,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_3.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.internal.gene.body.3 <- 
 apply(X = metageneplot.data.k4.internal.gene.body.3,
       MARGIN = 2,
       FUN = median,na.rm = T)

median.metageneplot.data.k4.internal.gene.body <- (median.metageneplot.data.k4.internal.gene.body.1 +
                                                    median.metageneplot.data.k4.internal.gene.body.2 +
                                                    median.metageneplot.data.k4.internal.gene.body.3)/3

plot(smooth.spline(x=1:40,y=median.metageneplot.data.k4.internal.gene.body,spar = 0.5),
     type="l",col="cyan",lwd=3,ylim=c(0,10),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)




metageneplot.data.k4.tss.1 <- 
 metageneplot(gene.names = k4.tss,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_1.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.tss.1 <- 
 apply(X = metageneplot.data.k4.tss.1,
       MARGIN = 2,
       FUN = median,na.rm = T)

metageneplot.data.k4.tss.2 <- 
 metageneplot(gene.names = k4.tss,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_2.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.tss.2 <- 
 apply(X = metageneplot.data.k4.tss.2,
       MARGIN = 2,
       FUN = median,na.rm = T)



metageneplot.data.k4.tss.3 <- 
 metageneplot(gene.names = k4.tss,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_3.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.tss.3 <- 
 apply(X = metageneplot.data.k4.tss.3,
       MARGIN = 2,
       FUN = median,na.rm = T)

median.metageneplot.data.k4.tss <- 
 (median.metageneplot.data.k4.tss.1 + 
   median.metageneplot.data.k4.tss.2 +
   median.metageneplot.data.k4.tss.3)/3

plot(smooth.spline(x=1:40,y=median.metageneplot.data.k4.tss,spar = 0.5),
     type="l",col="green",lwd=3,ylim=c(0,70),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)






metageneplot.data.k4.tes.1 <- 
 metageneplot(gene.names = k4.tes,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_1.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.tes.1 <- 
 apply(X = metageneplot.data.k4.tes.1,
       MARGIN = 2,
       FUN = median,na.rm = T)



metageneplot.data.k4.tes.2 <- 
 metageneplot(gene.names = k4.tes,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_2.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.tes.2 <- 
 apply(X = metageneplot.data.k4.tes.2,
       MARGIN = 2,
       FUN = median,na.rm = T)


metageneplot.data.k4.tes.3 <- 
 metageneplot(gene.names = k4.tes,
              genes = ostta.genes,
              bw.file.name = "chipseq_data/bw_data/h3k4me3/temp20/h3k4me3_20_chip_3.bw",
              ntile.txs = 10,
              ntile.gene.body = 20)
median.metageneplot.data.k4.tes.3 <- 
 apply(X = metageneplot.data.k4.tes.3,
       MARGIN = 2,
       FUN = median,na.rm = T)


median.metageneplot.data.k4.tes <- 
 (median.metageneplot.data.k4.tes.1 + 
   median.metageneplot.data.k4.tes.2 +
   median.metageneplot.data.k4.tes.3)/3

plot(smooth.spline(x=1:40,y=median.metageneplot.data.k4.tes,spar = 0.5),
     type="l",col="red",lwd=3,ylim=c(0,20),axes=F,xlab="",ylab="")
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)


jpeg(filename = "images/h3k4me3_metageneplots.jpg",res = 250,width = 900,height = 900)
plot(smooth.spline(x=1:40,y=median.metageneplot.data.k4.gene.body,spar = 0.5),
     type="l",col="blue",lwd=3,ylim=c(0,80),axes=F,xlab="",ylab="")
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k4.internal.gene.body,spar = 0.5),
      type="l",col="cyan",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k4.tss,spar = 0.5),
      type="l",col="green",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k4.tes,spar = 0.5),
      type="l",col="red",lwd=3)
axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)
dev.off()


jpeg(filename = "images/h3k24me3_expression_barplot_different_genes.jpg",res = 250,width = 600,height = 1000)
boxplot(gene.expression.20[k4.gene.body],
        gene.expression.20[k4.tss],
        gene.expression.20[k4.internal.gene.body],
        gene.expression.20[k4.tes],
        gene.expression.20[non.k4.genes],
        col=c("blue","green","cyan","red","grey"),
        outline=F)
dev.off()




library(clusterProfiler)
library(enrichplot)
library(org.Otauriv5.eg.db)
enrichment <- enrichGO(gene = all.k4.genes, pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                       OrgDb = org.Otauriv5.eg.db,
                       ont = "BP", 
                       keyType = "GID")
enrichment.df <- as.data.frame(enrichment)

write.table(x = enrichment.df, file = "go_enrichment_k4.tsv",quote = F,sep = "\t",row.names = F)

enrichment.df[which.max(enrichment.df$Count),]

sort(enrichment.df$Count,decreasing = T)[10:20]

enrichment.df[enrichment.df$Count == 26,]
enrichment.df[enrichment.df$Count == 24,]
enrichment.df[enrichment.df$Count == 22,]
enrichment.df[enrichment.df$Count == 21,]
enrichment.df[enrichment.df$Count == 20,]
enrichment.df[enrichment.df$Count == 18,]

enrichment.df[enrichment.df$Description == "positive regulation of biological process",]
enrichment.df[enrichment.df$Description == "locomotion",]


enrichment.df[enrichment.df$Description == "regulation of gene expression",]


barplot(c(gene.expression.20["ostta14g00070"],gene.expression.10["ostta14g00070"]))
barplot(c(gene.expression.20["ostta06g04460"],gene.expression.10["ostta06g04460"]))

gene.id <- "ostta11g02490"
barplot(c(gene.expression.20[gene.id],gene.expression.10[gene.id]))

jpeg(filename = "images/h3k27me3_treeplot.jpg",res = 250,width = 3600,height = 2000)
treeplot(x = pairwise_termsim(enrichment))
dev.off()

jpeg(filename = "../images/h3k27me3_treeplot_gene_body.jpg",res = 250,width = 3800,height = 1500)
treeplot(x = pairwise_termsim(k27.go.terms.filtered))
dev.off()

jpeg(filename = "images/h3k27me3_emapplot.jpg",res = 250,width = 2000,height = 2000)
emapplot(pairwise_termsim(enrichment))
dev.off()


cnetplot(pairwise_termsim(enrichment))

barplot(enrichment)

gene.ids <- ostta.genes$gene_id
k27.enrich.kegg <- enrichKEGG(gene  = paste0("OT_",all.k27.genes),
                                    universe = paste0("OT_",gene.ids),
                                    organism = "ota",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05)

df.activated.enrich.kegg <- as.data.frame(k27.enrich.kegg)
head(df.activated.enrich.kegg)

library(pathview)

pathview.genes <- rep(0,length(gene.ids))
names(pathview.genes) <- gene.ids

pathview.genes[all.k27.genes] <- 1

names(pathview.genes) <- paste0("OT_",gene.ids)

pathview(gene.data = sort(pathview.genes,decreasing = TRUE),
         pathway.id = "ota04814",
         species = "ota",
         gene.idtype = "KEGG")

library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006412","translation",4.38869169324396,41.635159951405136,0.6875841543305287,0,"translation"),
                     c("GO:0002181","cytoplasmic translation",0.3974763229665376,8.595948611417302,0.755633504893879,0.6792736,"translation"),
                     c("GO:0009059","macromolecule biosynthetic process",16.15219652221219,32.13388504723149,0.8406622568296607,0.37915401,"translation"),
                     c("GO:0043603","amide metabolic process",6.707376287050344,36.45011550587449,0.9109747285969385,0.18706276,"translation"),
                     c("GO:0009765","photosynthesis, light harvesting",0.02241434707504848,4.839062538560749,0.952328578586525,0.07239153,"photosynthesis, light harvesting"),
                     c("GO:0015979","photosynthesis",0.228607115192195,7.075542911759376,0.941812778414762,0.0906509,"photosynthesis"),
                     c("GO:0048518","positive regulation of biological process",4.759580864033397,2.926034087880324,0.8707557167557012,-0,"positive regulation of biological process"),
                     c("GO:0010557","positive regulation of macromolecule biosynthetic process",2.093180586483937,2.989915809683905,0.7484373322239899,0.27791144,"positive regulation of biological process"),
                     c("GO:0031323","regulation of cellular metabolic process",14.819661863282441,3.0711442275048717,0.7768013068864807,0.54794604,"positive regulation of biological process"),
                     c("GO:0140694","non-membrane-bounded organelle assembly",0.7232459602662196,6.425049931933804,0.6134007334542475,0.01308302,"ribosome biogenesis"),
                     c("GO:0000027","ribosomal large subunit assembly",0.012821578346646382,3.275341540820152,0.6307593041874131,0.61078702,"ribosome biogenesis"),
                     c("GO:0022613","ribonucleoprotein complex biogenesis",2.5375929564430133,4.856610166214762,0.6332784349267674,0.65931649,"ribosome biogenesis"),
                     c("GO:0042254","ribosome biogenesis",2.136224808753416,5.364416283328431,0.566141578351784,0.5650077,"ribosome biogenesis"),
                     c("GO:0070925","organelle assembly",1.3561320812847646,3.9914667506842227,0.6213427911227302,0.65542123,"ribosome biogenesis"),
                     c("GO:0071826","protein-RNA complex organization",0.6518695497816881,3.3975782710402287,0.7204190675227414,0.36838395,"ribosome biogenesis"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
jpeg(filename = "images/h3k4me3_revigo_treemap.jpg",res = 300,width = 1400,height = 900)

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


overlap.k4.k27 <- read.table(file="chipseq_data/bed_data/h3k27me3_h3k4me3_ll_20.bed",header=F,sep = "\t",as.is = T)

max(overlap.k4.k27$V3 - overlap.k4.k27$V2)
summary(overlap.k4.k27$V3 - overlap.k4.k27$V2)
nrow(overlap.k4.k27)


k4.k27.genes <- intersect(all.k27.genes, all.k4.genes)
summary(gene.expression.20[k4.k27.genes])
summary(gene.expression.20[all.k4.genes])
summary(gene.expression.20[all.k27.genes])


jpeg(filename = "images/bivalent_expression_barplot.jpg",res = 250,width = 600,height = 1000)
boxplot(gene.expression.20[all.k27.genes],
        gene.expression.20[non.k4.k27.genes],
        gene.expression.20[k4.k27.genes],
        gene.expression.20[all.k4.genes],
        col=c("blue","grey",colorRampPalette(c('blue','darkorange'))(3)[2],"darkorange"),outline=F)
dev.off()
wilcox.test(x = gene.expression.20[k4.k27.genes], y = gene.expression.20[non.k4.k27.genes])

non.k4.k27.genes <- intersect(non.k27.genes,non.k4.genes)



sum(gene.expression.20[k4.k27.genes] > 1)/length(k4.k27.genes)
sum(gene.expression.20[k4.k27.genes] > 10)/length(k4.k27.genes)

sum(gene.expression.20[k4.k27.genes] < 10)/length(k4.k27.genes)
sum(gene.expression.20[k4.k27.genes] >= 10)/length(k4.k27.genes)

sum(gene.expression.20[all.k27.genes] > 10)/length(all.k27.genes)


library(VennDiagram)

jpeg(filename = "images/h3k27me3_h3k4me3_venndiagram.jpg",res = 250,width = 1500,height = 1500)
grid.newpage()
draw.pairwise.venn(area1 = length(all.k27.genes),cat.pos = c(180,180),cat.dist = 0.05,
                   area2 = length(all.k4.genes),
                   cross.area = length(intersect(all.k27.genes,all.k4.genes)),
                   lwd = 3,category = c("H3K27me3","H3K4me3"),euler.d = T,
                   col = c("blue","darkorange"),
                   fill = c(" blue","darkorange"),alpha = 0.5,
                   cex = 2,
                   cat.cex = 2)
dev.off()


head(ll.20.peak.annotation)

table(subset(ll.20.peak.annotation, target.genes %in% k4.k27.genes)[["peak.annotation"]])

k4.k27.gene.annotation <- subset(ll.20.peak.annotation, target.genes %in% k4.k27.genes)
100*table(k4.k27.gene.annotation[!duplicated(k4.k27.gene.annotation$target.genes),][["peak.annotation"]])/nrow(k4.k27.gene.annotation)

enrichment <- enrichGO(gene = k4.k27.genes, pvalueCutoff = 1,qvalueCutoff = 1,
                       OrgDb = org.Otauriv5.eg.db,
                       ont = "BP", 
                       keyType = "GID")
enrichment.df <- as.data.frame(enrichment)

max(enrichment.df$Count)
sort(enrichment.df$Count)

which.max(enrichment.df$Count)
enrichment.df[3,]
enrichment.df[enrichment.df$Count == 3,]




gene.expression.20 <- 
 (temp.gene.expression$T20_1 + 
   temp.gene.expression$T20_2 + 
   temp.gene.expression$T20_3)/3


genes.1 <- c("ostta09g03765","ostta02g05090","ostta04g03540","ostta12g02370")
genes.2 <- c("ostta05g01050","ostta16g01360")

genes.k27.k4.expression.1 <- temp.gene.expression[genes.1,c("T20_1","T20_2","T20_3")]
genes.k27.k4.expression.2 <- temp.gene.expression[genes.2,c("T20_1","T20_2","T20_3")]


genes.mean.1 <- apply(X = genes.k27.k4.expression.1,MARGIN = 1,FUN = mean)
genes.sd.1 <- apply(X = genes.k27.k4.expression.1,MARGIN = 1,FUN = sd)

genes.mean.2 <- apply(X = genes.k27.k4.expression.2,MARGIN = 1,FUN = mean)
genes.sd.2 <- apply(X = genes.k27.k4.expression.2,MARGIN = 1,FUN = sd)

jpeg(filename = "images/bivalent_example_repressed.jpg",res = 250,width = 1500,height = 1500)
par(mar=c(8,5,2,2),lwd=1.5)
xpos <- barplot(genes.mean.1,col=colorRampPalette(c('blue','darkorange'))(3)[2],ylim=c(0,12),las=2)
arrows(x0 = xpos,y0 = genes.mean.1+genes.sd.1,x1 = xpos, y1 = genes.mean.1-genes.sd.1,code=3,angle=90,length=0.05)
dev.off()

jpeg(filename = "images/bivalent_example_activated.jpg",res = 250,width = 750,height = 1500)
par(mar=c(8,5,2,2),lwd=1.5)
barplot(genes.mean.2,col=colorRampPalette(c('blue','darkorange'))(3)[2],las=2,ylim=c(0,600))
arrows(x0 = xpos,y0 = genes.mean.2+genes.sd.2,x1 = xpos, y1 = genes.mean.2-genes.sd.2,code=3,angle=90,length=0.05)
dev.off()


library(rtracklayer)
file.1 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw"
file.2 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_2n.bw"
file.3 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_3n.bw"

file.1 <- "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_1.bw"
file.2 <- "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_2.bw"
file.3 <- "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_1.bw"


bw.file.1 <- import.bw(con = file.1)
bw.data.1 <- as.data.frame(bw.file.1)

bw.file.2 <- import.bw(con = file.2)
bw.data.2 <- as.data.frame(bw.file.2)

bw.file.3 <- import.bw(con = file.3)
bw.data.3 <- as.data.frame(bw.file.3)

bw.data.1.1 <- subset(bw.data.1, seqnames==16)
bw.data.1.2 <- subset(bw.data.2, seqnames==16)
bw.data.1.3 <- subset(bw.data.3, seqnames==16)
signal <- vector(mode = "numeric", length = floor(chr.lens[16]/200))

ini <- 1
stop <- ini + 200
j <- 1
while(stop < chr.lens[16])
{
 print(j)
 signal[j] <- (mean(subset(x = bw.data.1.1, 
                                 start >= ini & 
                                  end < stop)$score) +
                       mean(subset(x = bw.data.1.2, 
                                   start >= ini & 
                                    end < stop)$score) +
                       mean(subset(x = bw.data.1.3, 
                                   start >= ini & 
                                    end < stop)$score))/3
 j <- j + 1
 ini <- stop
 stop <- ini + 200
}


plot(signal,type="l")

signal.20 <- signal
signal.zt8 <- signal

plot(signal.20,type="l")
lines(signal.zt8,type="l",col="lightgrey")


cca1.zt0.bw <- as.data.frame(import.bw(con = "/home/fran/Nextcloud2/Microalgas/electra/ostreococcus/chip_rnaseq_cca1_sd_20/CCA1/SD_20oC_ZT0_1_IP.bw"))
cca1.zt4.bw <- as.data.frame(import.bw(con = "/home/fran/Nextcloud2/Microalgas/electra/ostreococcus/chip_rnaseq_cca1_sd_20/CCA1/SD_20oC_ZT4_1_IP.bw"))
cca1.zt8.bw <- as.data.frame(import.bw(con = "/home/fran/Nextcloud2/Microalgas/electra/ostreococcus/chip_rnaseq_cca1_sd_20/CCA1/SD_20oC_ZT8_1_IP.bw"))
cca1.zt12.bw <- as.data.frame(import.bw(con = "/home/fran/Nextcloud2/Microalgas/electra/ostreococcus/chip_rnaseq_cca1_sd_20/CCA1/SD_20oC_ZT12_1_IP.bw"))
cca1.zt16.bw <- as.data.frame(import.bw(con = "/home/fran/Nextcloud2/Microalgas/electra/ostreococcus/chip_rnaseq_cca1_sd_20/CCA1/SD_20oC_ZT16_1_IP.bw"))
cca1.zt20.bw <- as.data.frame(import.bw(con = "/home/fran/Nextcloud2/Microalgas/electra/ostreococcus/chip_rnaseq_cca1_sd_20/CCA1/SD_20oC_ZT20_1_IP.bw"))

chr <- 6
ini <- 198264
stop <- 199112

current.peak.zt0 <- max(subset(cca1.zt0.bw, seqnames == chr & start >= ini & end <= stop)$score)
current.peak.zt4 <- max(subset(cca1.zt4.bw, seqnames == chr & start >= ini & end <= stop)$score)
current.peak.zt8 <- max(subset(cca1.zt8.bw, seqnames == chr & start >= ini & end <= stop)$score)
current.peak.zt12 <- max(subset(cca1.zt12.bw, seqnames == chr & start >= ini & end <= stop)$score)
current.peak.zt16 <- max(subset(cca1.zt16.bw, seqnames == chr & start >= ini & end <= stop)$score)
current.peak.zt20 <- max(subset(cca1.zt20.bw, seqnames == chr & start >= ini & end <= stop)$score)

barplot(height = c(current.peak.zt0,current.peak.zt4,current.peak.zt8,current.peak.zt12, current.peak.zt16, current.peak.zt20))

sd.20.gene.expression <- read.table(file="/home/fran/Nextcloud2/Microalgas/electra/ostreococcus/chip_rnaseq_cca1_sd_20/CCA1/sd_20_gene_expression.txt",header=T,as.is=T)
head(sd.20.gene.expression)
current.gene <- "ostta06g02940"
current.gene <- "ostta07g03440"
current.gene < "ostta18g01570"
current.gene <- "ostta13g01820"
current.gene <- "ostta06g01220"

current.gene.expression <- (unlist(sd.20.gene.expression[current.gene,])[1:6] +
unlist(sd.20.gene.expression[current.gene,])[7:12] +
unlist(sd.20.gene.expression[current.gene,])[13:18])/3

plot(current.gene.expression,type="l")
lines(current.gene.expression,type="l")






# Entering data 
zt <- c(0, 4, 8, 12, 16, 20) 
cca1.binding <- c(current.peak.zt0,current.peak.zt4,current.peak.zt8,current.peak.zt12, current.peak.zt16, current.peak.zt20) 
current.gene.expression#c(0.3, 0.25, 0.3, 0.5, 0.4, 0.2, 0.6) 

# Creating Data Frame 
perf <- data.frame(zt, cca1.binding, current.gene.expression) 

# Plotting Charts and adding a secondary axis 
library(ggplot2) 
ggp <- ggplot(perf)  +  
 geom_bar(aes(x=zt, y=cca1.binding),stat="identity", fill="cyan",colour="#006000")+ 
 geom_line(aes(x=zt, y=3*current.gene.expression),stat="identity",color="red",size=2)+ 
 labs(title= "CCA1 binding vs Gene Expression", 
      x="Time (ZT)",y="CCA1 Binding Signal (CPM)")+ 
 scale_y_continuous(sec.axis=sec_axis(~./3,name="Gene Expression (FPKM)")) 
ggp 

## 

mapped.reads.k27.20c.1 <- 3146662/1E6 
mapped.reads.k27.20c.2 <- 5697849/1E6
mapped.reads.k27.20c.3 <- 5157582/1E6

mapped.reads.k27.20c.zt8.1 <- 2254676/1E6 
mapped.reads.k27.20c.zt8.2 <- 2639017/1E6
mapped.reads.k27.20c.zt16.1 <- 3605883/1E6
mapped.reads.k27.20c.zt16.2 <- 3405826/1E6


read.counts.k27.temp20c.1 <- read.table(file ="chipseq_data/consensus_counts/h3k27me3_ll_ld_consensus_ll1_counts.bed", header = F, sep = "\t" )
head(read.counts.k27.temp20c.1)
nrow(read.counts.k27.temp20c.1)
cpm.k27.20.1 <- read.counts.k27.temp20c.1$V10 / mapped.reads.k27.20c.1

read.counts.k27.temp20c.2 <- read.table(file = "chipseq_data/consensus_counts/h3k27me3_ll_ld_consensus_ll2_counts.bed", header = F, sep = "\t" )
head(read.counts.k27.temp20c.2)
nrow(read.counts.k27.temp20c.2)
cpm.k27.20.2 <- read.counts.k27.temp20c.2$V10 / mapped.reads.k27.20c.2

read.counts.k27.temp20c.3 <- read.table(file = "chipseq_data/consensus_counts/h3k27me3_ll_ld_consensus_ll3_counts.bed", header = F, sep = "\t" )
head(read.counts.k27.temp20c.3)
nrow(read.counts.k27.temp20c.3)
cpm.k27.20.3 <- read.counts.k27.temp20c.3$V10 / mapped.reads.k27.20c.3

read.counts.k27.20c.zt8.1 <- read.table(file = "chipseq_data/consensus_counts/h3k27me3_ll_ld_consensus_zt8_1_counts.bed", header = F, sep = "\t" )
head(read.counts.k27.20c.zt8.1)
nrow(read.counts.k27.20c.zt8.1)
cpm.k27.zt8.1 <- read.counts.k27.20c.zt8.1$V10 / mapped.reads.k27.20c.zt8.1

read.counts.k27.20c.zt8.2 <- read.table(file = "chipseq_data/consensus_counts/h3k27me3_ll_ld_consensus_zt8_2_counts.bed", header = F, sep = "\t" )
head(read.counts.k27.20c.zt8.2)
nrow(read.counts.k27.20c.zt8.2)
cpm.k27.zt8.2 <- read.counts.k27.20c.zt8.2$V10 / mapped.reads.k27.20c.zt8.2

read.counts.k27.20c.zt16.1 <- read.table(file = "chipseq_data/consensus_counts/h3k27me3_ll_ld_consensus_zt16_1_counts.bed", header = F, sep = "\t" )
head(read.counts.k27.20c.zt16.1)
nrow(read.counts.k27.20c.zt16.1)
cpm.k27.zt16.1 <- read.counts.k27.20c.zt16.1$V10 / mapped.reads.k27.20c.zt16.1

read.counts.k27.20c.zt16.2 <- read.table(file = "chipseq_data/consensus_counts/h3k27me3_ll_ld_consensus_zt16_2_counts.bed", header = F, sep = "\t" )
head(read.counts.k27.20c.zt16.2)
nrow(read.counts.k27.20c.zt16.2)
cpm.k27.zt16.2 <- read.counts.k27.20c.zt16.2$V10 / mapped.reads.k27.20c.zt16.2

head(read.counts.k27.20c.zt16.1[,1:3])
peak.names <- apply(X = read.counts.k27.20c.zt16.1[,1:3],MARGIN = 1,FUN = paste,collapse="_")

cpm.data.ll.zt8.zt16 <- matrix(data = c(cpm.k27.20.1,cpm.k27.20.2,cpm.k27.20.3,cpm.k27.zt8.1,cpm.k27.zt8.2,cpm.k27.zt16.1,cpm.k27.zt16.2),ncol=7)
rownames(cpm.data.ll.zt8.zt16) <- peak.names
colnames(cpm.data.ll.zt8.zt16) <- c("LL20_1","LL20_2","LL20_3","LDZT8_1","LDZT8_2","LDZT16_1","LDZT16_2")
boxplot(cpm.data.ll.zt8.zt16,outline=F,las=2,
        col=c(rep("blue",3),rep("yellow",2),rep("grey",2)))

write.table(x = cpm.data.ll.zt8.zt16, file = "cpm_data_ll_zt8_zt16.tsv", 
            sep = "\t", quote = F,row.names = F)

library(NormalyzerDE)

design <- data.frame(sample=colnames(cpm.data.ll.zt8.zt16),
                     group=c(rep("LL20",3),rep("LDZT8",2),rep("LDZT16",2)))

write.table(x = design,file = "normalyzer_design_ll20_ldzt8_ldzt16.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "ll20_ldzt8_ldzt16",designPath = "normalyzer_design_ll20_ldzt8_ldzt16.tsv",
           dataPath = "cpm_data_ll_zt8_zt16.tsv",outputDir = ".")


normalized.log2.cpm <- read.table(file="ll20_ldzt8_ldzt16/Quantile-normalized.txt", header=T)
rownames(normalized.log2.cpm) <- peak.names
head(normalized.log2.cpm)

boxplot(normalized.log2.cpm,outline=F,las=2,
        col=c(rep("blue",3),rep("yellow",2),rep("grey",2)))


plot(x = normalized.log2.cpm[,1],normalized.log2.cpm[,2])
plot(x = normalized.log2.cpm[,1],normalized.log2.cpm[,3])
plot(x = normalized.log2.cpm[,2],normalized.log2.cpm[,3])

plot(x = normalized.log2.cpm[,4],normalized.log2.cpm[,5])

plot(x = normalized.log2.cpm[,6],normalized.log2.cpm[,7])

ll.20 <- rowMeans(normalized.log2.cpm[,1:3])
zt8 <- rowMeans(normalized.log2.cpm[,4:5])
zt16 <- rowMeans(normalized.log2.cpm[,6:7])

plot(ll.20,zt8,pch=19)
lines(x=c(0,20),y=c(0,20),col="red",lwd=2)

plot(ll.20,zt16,pch=19)
lines(x=c(0,20),y=c(0,20),col="red",lwd=2)

plot(zt8,zt16,pch=19)
lines(x=c(0,20),y=c(0,20),col="red",lwd=2)

which((ll.20 - zt8) > 1)
which((ll.20 - zt8) < -1)

which((ll.20 - zt16) > 1)
which((ll.20 - zt16) < -1)

length(which((ll.20 - zt8) > 1))
length(which((ll.20 - zt16) > 1))
length(intersect(which((ll.20 - zt8) > 1),which((ll.20 - zt16) > 1)))

library(limma)
limma.experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,3,3)))
colnames(limma.experimental.design) <- c("LL20", "LDZT8", "LDZT16")

linear.fit <- lmFit(normalized.log2.cpm, limma.experimental.design)

contrast.matrix <- makeContrasts(LL20-LDZT8, LL20-LDZT16, LDZT16-LDZT8, levels = c("LL20", "LDZT8", "LDZT16"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)                               

ll20.zt8 <- topTable(contrast.results, 
                     number=nrow(normalized.log2.cpm), 
                     coef = 1, sort.by = "logFC", )
head(ll20.zt8)

ll20.zt8.logfc <- ll20.zt8$logFC
ll20.zt8.q.val <- ll20.zt8$adj.P.Val

activated.ll20.zt8 <- rownames(ll20.zt8)[ll20.zt8.logfc > 1 & ll20.zt8.q.val < 0.05]
length(activated.ll20.zt8)

repressed.ll20.zt8 <- rownames(ll20.zt8)[ll20.zt8.logfc < -1 & ll20.zt8.q.val < 0.05]
length(repressed.ll20.zt8)


ll20.zt16 <- topTable(contrast.results, 
                      number=nrow(normalized.log2.cpm), 
                      coef = 2, sort.by = "logFC", )
head(ll20.zt16)

ll20.zt16.logfc <- ll20.zt16$logFC
ll20.zt16.q.val <- ll20.zt16$adj.P.Val

activated.ll20.zt16 <- rownames(ll20.zt16)[ll20.zt16.logfc > 1 & ll20.zt16.q.val < 0.05]
length(activated.ll20.zt16)

repressed.ll20.zt16 <- rownames(ll20.zt16)[ll20.zt16.logfc < -1 & ll20.zt16.q.val < 0.05]
length(repressed.ll20.zt16)

length(intersect(activated.ll20.zt8,activated.ll20.zt16))
2_161795_165554
ostta10g00130
length(intersect(repressed.ll20.zt8,repressed.ll20.zt16))

temp.cols <- c("#00FFFF","#6699FF","#CC33FF","#FF0000")
names(temp.cols) <- c("10C","14C","20C","26C")
#col=c("#00FFFF","#6699FF","#CC33FF","#FF0000")
#"10ºC","14ºC","20ºC","26ºC"
#11_333391_333724 ostta11g01770
#This family of proteins includes secreted effectors with a role in host-pathogen interactions. One member is known to contain an RxLR-dEER motif, which is involved in the translocation of the effector into the host cell by binding to specific phospholipids on the cell surface. The same member is also implicated in the suppression of host cell death triggered by other proteins, indicating a possible function in manipulating host cell pathways to benefit the pathogen. Other members are annotated as putative ankyrin repeat proteins, suggesting a potential role in protein-protein interactions, although their specific functions remain uncharacterized.


mapped.reads.k4.20c.1 <- 6832850/1E6 
mapped.reads.k4.20c.2 <- 5993486/1E6
mapped.reads.k4.20c.3 <- 6701422/1E6
mapped.reads.k4.20c.zt8.1 <- 2191204/1E6 
mapped.reads.k4.20c.zt8.2 <- 3898985/1E6
mapped.reads.k4.20c.zt16.1 <- 4894407/1E6
mapped.reads.k4.20c.zt16.2 <- 3797062/1E6

read.counts.k4.temp20c.1 <- read.table(file ="chipseq_data/consensus_counts/h3k4me3_ll_ld_consensus_ll1_counts.bed", header = F, sep = "\t" )
head(read.counts.k4.temp20c.1)
nrow(read.counts.k4.temp20c.1)
cpm.k4.20.1 <- read.counts.k4.temp20c.1$V11 / mapped.reads.k4.20c.1

read.counts.k4.temp20c.2 <- read.table(file ="chipseq_data/consensus_counts/h3k4me3_ll_ld_consensus_ll2_counts.bed", header = F, sep = "\t" )
head(read.counts.k4.temp20c.2)
nrow(read.counts.k4.temp20c.2)
cpm.k4.20.2 <- read.counts.k4.temp20c.2$V11 / mapped.reads.k4.20c.2

read.counts.k4.temp20c.3 <- read.table(file ="chipseq_data/consensus_counts/h3k4me3_ll_ld_consensus_ll3_counts.bed", header = F, sep = "\t" )
head(read.counts.k4.temp20c.3)
nrow(read.counts.k4.temp20c.3)
cpm.k4.20.3 <- read.counts.k4.temp20c.3$V11 / mapped.reads.k4.20c.3

read.counts.k4.20c.zt8.1 <- read.table(file = "chipseq_data/consensus_counts/h3k4me3_ll_ld_consensus_zt8_1_counts.bed", header = F, sep = "\t" )
head(read.counts.k4.20c.zt8.1)
nrow(read.counts.k4.20c.zt8.1)
cpm.k4.zt8.1 <- read.counts.k4.20c.zt8.1$V11 / mapped.reads.k4.20c.zt8.1

read.counts.k4.20c.zt8.2 <- read.table(file = "chipseq_data/consensus_counts/h3k4me3_ll_ld_consensus_zt8_2_counts.bed", header = F, sep = "\t" )
head(read.counts.k4.20c.zt8.2)
nrow(read.counts.k4.20c.zt8.2)
cpm.k4.zt8.2 <- read.counts.k4.20c.zt8.2$V11 / mapped.reads.k4.20c.zt8.2

read.counts.k4.20c.zt16.1 <- read.table(file = "chipseq_data/consensus_counts/h3k4me3_ll_ld_consensus_zt16_1_counts.bed", header = F, sep = "\t" )
head(read.counts.k4.20c.zt16.1)
nrow(read.counts.k4.20c.zt16.1)
cpm.k4.zt16.1 <- read.counts.k4.20c.zt16.1$V11 / mapped.reads.k4.20c.zt16.1

read.counts.k4.20c.zt16.2 <- read.table(file = "chipseq_data/consensus_counts/h3k4me3_ll_ld_consensus_zt16_2_counts.bed", header = F, sep = "\t" )
head(read.counts.k4.20c.zt16.2)
nrow(read.counts.k4.20c.zt16.2)
cpm.k4.zt16.2 <- read.counts.k4.20c.zt16.2$V11 / mapped.reads.k4.20c.zt16.2

k4.peak.names <- apply(X = read.counts.k4.20c.zt16.1[,1:3],MARGIN = 1,FUN = paste,collapse="_")

k4.cpm.data.ll.zt8.zt16 <- matrix(data = c(cpm.k4.20.1,
                                           cpm.k4.20.2,
                                           cpm.k4.20.3,
                                           cpm.k4.zt8.1,
                                           cpm.k4.zt8.2,
                                           cpm.k4.zt16.1,
                                           cpm.k4.zt16.2),ncol=7)
rownames(k4.cpm.data.ll.zt8.zt16) <- k4.peak.names
colnames(k4.cpm.data.ll.zt8.zt16) <- c("LL20_1","LL20_2","LL20_3","LDZT8_1","LDZT8_2","LDZT16_1","LDZT16_2")
boxplot(k4.cpm.data.ll.zt8.zt16,outline=F,las=2,
        col=c(rep(temp.cols["20C"],3),rep("yellow",2),rep("grey",2)))

write.table(x = k4.cpm.data.ll.zt8.zt16, file = "k4_cpm_data_ll_zt8_zt16.tsv", 
            sep = "\t", quote = F,row.names = F)

library(NormalyzerDE)

design <- data.frame(sample=colnames(k4.cpm.data.ll.zt8.zt16),
                     group=c(rep("LL20",3),rep("LDZT8",2),rep("LDZT16",2)))

write.table(x = design,file = "normalyzer_design_k4_ll20_ldzt8_ldzt16.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "k4_ll20_ldzt8_ldzt16",designPath = "normalyzer_design_k4_ll20_ldzt8_ldzt16.tsv",
           dataPath = "k4_cpm_data_ll_zt8_zt16.tsv",outputDir = ".")


k4.normalized.log2.cpm <- read.table(file="k4_ll20_ldzt8_ldzt16/Quantile-normalized.txt", header=T)
rownames(k4.normalized.log2.cpm) <- k4.peak.names
head(k4.normalized.log2.cpm)

boxplot(k4.normalized.log2.cpm,outline=F,las=2,
        col=c(rep(temp.cols["20C"],3),rep("yellow",2),rep("grey",2)))


plot(x = k4.normalized.log2.cpm[,1],k4.normalized.log2.cpm[,2])
plot(x = k4.normalized.log2.cpm[,1],k4.normalized.log2.cpm[,3])
plot(x = k4.normalized.log2.cpm[,2],k4.normalized.log2.cpm[,3])

plot(x = k4.normalized.log2.cpm[,4],k4.normalized.log2.cpm[,5])

plot(x = k4.normalized.log2.cpm[,6],k4.normalized.log2.cpm[,7])

plot(x = k4.normalized.log2.cpm[,4],k4.normalized.log2.cpm[,6])

k4.ll.20 <- rowMeans(k4.normalized.log2.cpm[,1:3])
k4.zt8 <- rowMeans(k4.normalized.log2.cpm[,4:5])
k4.zt16 <- rowMeans(k4.normalized.log2.cpm[,6:7])

plot(k4.ll.20,k4.zt8,pch=19)
lines(x=c(0,20),y=c(0,20),col="red",lwd=2)

plot(k4.ll.20,k4.zt16,pch=19)
lines(x=c(0,20),y=c(0,20),col="red",lwd=2)

plot(k4.zt8,k4.zt16,pch=19)
lines(x=c(0,20),y=c(0,20),col="red",lwd=2)

which((k4.ll.20 - k4.zt8) > 1)
which((k4.ll.20 - k4.zt8) < -1)

which((k4.ll.20 - k4.zt16) > 1)
which((ll.20 - zt16) < -1)

length(which((ll.20 - zt8) > 1))
length(which((ll.20 - zt16) > 1))
length(intersect(which((ll.20 - zt8) > 1),which((ll.20 - zt16) > 1)))

library(limma)
k4.limma.experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,3,3)))
colnames(k4.limma.experimental.design) <- c("LL20", "LDZT8", "LDZT16")

k4.linear.fit <- lmFit(k4.normalized.log2.cpm, k4.limma.experimental.design)

k4.contrast.matrix <- makeContrasts(LL20-LDZT8, LL20-LDZT16, LDZT16-LDZT8, levels = c("LL20", "LDZT8", "LDZT16"))

k4.contrast.linear.fit <- contrasts.fit(k4.linear.fit, k4.contrast.matrix)
k4.contrast.results <- eBayes(k4.contrast.linear.fit)                               

k4.ll20.zt8 <- topTable(k4.contrast.results, 
                     number=nrow(k4.normalized.log2.cpm), 
                     coef = 1, sort.by = "logFC", )
head(k4.ll20.zt8)

k4.ll20.zt8.logfc <- k4.ll20.zt8$logFC
k4.ll20.zt8.q.val <- k4.ll20.zt8$adj.P.Val

k4.activated.ll20.zt8 <- rownames(k4.ll20.zt8)[k4.ll20.zt8.logfc > 1 & k4.ll20.zt8.q.val < 0.05]
length(k4.activated.ll20.zt8)

k4.repressed.ll20.zt8 <- rownames(k4.ll20.zt8)[k4.ll20.zt8.logfc < -1 & k4.ll20.zt8.q.val < 0.05]
length(k4.repressed.ll20.zt8)


k4.ll20.zt16 <- topTable(k4.contrast.results, 
                      number=nrow(k4.normalized.log2.cpm), 
                      coef = 2, sort.by = "logFC", )
head(k4.ll20.zt16)

k4.ll20.zt16.logfc <- k4.ll20.zt16$logFC
k4.ll20.zt16.q.val <- k4.ll20.zt16$adj.P.Val

k4.activated.ll20.zt16 <- rownames(k4.ll20.zt16)[k4.ll20.zt16.logfc > 1 & k4.ll20.zt16.q.val < 0.05]
length(k4.activated.ll20.zt16)

k4.repressed.ll20.zt16 <- rownames(k4.ll20.zt16)[k4.ll20.zt16.logfc < -1 & k4.ll20.zt16.q.val < 0.05]
length(k4.repressed.ll20.zt16)

length(intersect(k4.activated.ll20.zt8,k4.activated.ll20.zt16))
2_161795_165554
ostta10g00130
length(intersect(k4.repressed.ll20.zt8,k4.repressed.ll20.zt16))


library(rtracklayer)
file.1 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw"
file.2 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_2n.bw"
file.3 <- "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_3n.bw"

file.1 <- "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_1.bw"
file.2 <- "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_2.bw"
file.3 <- "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_1.bw"


signal.chr <- function(file.1,file.2,file.3,chr.lens,chr,interval)
{
 bw.file.1 <- import.bw(con = file.1)
 bw.data.1 <- as.data.frame(bw.file.1)
 
 bw.file.2 <- import.bw(con = file.2)
 bw.data.2 <- as.data.frame(bw.file.2)
 
 bw.file.3 <- import.bw(con = file.3)
 bw.data.3 <- as.data.frame(bw.file.3)
 
 bw.data.1.1 <- subset(bw.data.1, seqnames==16)
 bw.data.1.2 <- subset(bw.data.2, seqnames==16)
 bw.data.1.3 <- subset(bw.data.3, seqnames==16)
 signal <- vector(mode = "numeric", length = floor(chr.lens[chr]/interval))
 
 ini <- 1
 stop <- ini + interval
 j <- 1
 while(stop < chr.lens[chr])
 {
  print(j)
  signal[j] <- (mean(subset(x = bw.data.1.1, 
                            start >= ini & 
                             end < stop)$score) +
                 mean(subset(x = bw.data.1.2, 
                             start >= ini & 
                              end < stop)$score) +
                 mean(subset(x = bw.data.1.3, 
                             start >= ini & 
                              end < stop)$score))/3
  j <- j + 1
  ini <- stop
  stop <- ini + interval
 }
 
 return(signal)
}


file.1 = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",
file.2 = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_2n.bw",
file.3 = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_3n.bw"


signal.k27.20 <- signal.chr(file.1 = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",
           file.2 = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",
           file.3 = "chipseq_data/bw_data/h3k27me3/temp20/h3k27me3_20_chip_1.bw",
           chr.lens = chr.lens,chr = 10,interval = 200)

signal.k27.zt8 <- signal.chr(file.1 = "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_1.bw",
                             file.2 = "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_1.bw",
                             file.3 = "chipseq_data/bw_data/h3k27me3/ld_20_zt8/h3k27me3_20_zt8_1.bw",
                             chr.lens = chr.lens,chr = 10,interval = 200)

signal.k27.zt16 <- signal.chr(file.1 = "chipseq_data/bw_data/h3k27me3/ld_20_zt16/h3k27me3_20_zt16_1.bw",
                             file.2 = "chipseq_data/bw_data/h3k27me3/ld_20_zt16/h3k27me3_20_zt16_1.bw",
                             file.3 = "chipseq_data/bw_data/h3k27me3/ld_20_zt16/h3k27me3_20_zt16_1.bw",
                             chr.lens = chr.lens,chr = 10,interval = 200)

bw.file.1 <- import.bw(con = file.1)
bw.data.1 <- as.data.frame(bw.file.1)

bw.file.2 <- import.bw(con = file.2)
bw.data.2 <- as.data.frame(bw.file.2)

bw.file.3 <- import.bw(con = file.3)
bw.data.3 <- as.data.frame(bw.file.3)

bw.data.1.1 <- subset(bw.data.1, seqnames==16)
bw.data.1.2 <- subset(bw.data.2, seqnames==16)
bw.data.1.3 <- subset(bw.data.3, seqnames==16)
signal <- vector(mode = "numeric", length = floor(chr.lens[16]/200))

ini <- 1
stop <- ini + 200
j <- 1
while(stop < chr.lens[16])
{
 print(j)
 signal[j] <- (mean(subset(x = bw.data.1.1, 
                           start >= ini & 
                            end < stop)$score) +
                mean(subset(x = bw.data.1.2, 
                            start >= ini & 
                             end < stop)$score) +
                mean(subset(x = bw.data.1.3, 
                            start >= ini & 
                             end < stop)$score))/3
 j <- j + 1
 ini <- stop
 stop <- ini + 200
}


plot(signal,type="l")

signal.20 <- signal
signal.zt8 <- signal

plot(signal.k27.20,type="l",col=temp.cols["20C"])
lines(signal.k27.zt8,type="l",col="yellow")
lines(signal.k27.zt8,type="l",col="grey")









chip.signal.zt8.1 <- read.table(
  file = "chip_data/k27_consensus_ld_20/h3k27me3_ld_20_zt8_1_counts.bed",
  header = F,sep = "\t")
chip.signal.zt8.1[,4] <- chip.signal.zt8.1[,4]/mapped.reads.k27.20c.zt8.1

chip.signal.zt8.2 <- read.table(
  file = "chip_data/k27_consensus_ld_20/h3k27me3_ld_20_zt8_2_counts.bed",
  header = F,sep = "\t")
chip.signal.zt8.2[,4] <- chip.signal.zt8.2[,4]/mapped.reads.k27.20c.zt8.2

plot(log2(chip.signal.zt8.1[,4]),log2(chip.signal.zt8.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

chip.signal.zt16.1 <- read.table(
  file = "chip_data/k27_consensus_ld_20/h3k27me3_ld_20_zt16_1_counts.bed",
  header = F,sep = "\t")
chip.signal.zt16.1[,4] <- chip.signal.zt16.1[,4]/mapped.reads.k27.20c.zt16.1

chip.signal.zt16.2 <- read.table(
  file = "chip_data/k27_consensus_ld_20/h3k27me3_ld_20_zt16_2_counts.bed",
  header = F,sep = "\t")
chip.signal.zt16.2[,4] <- chip.signal.zt16.2[,4]/mapped.reads.k27.20c.zt16.2

plot(log2(chip.signal.zt16.1[,4]),log2(chip.signal.zt16.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

plot(log2(chip.signal.zt8.1[,4]),log2(chip.signal.zt16.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

plot(log2(chip.signal.zt8.1[,4]),log2(chip.signal.zt16.1[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

plot(log2(chip.signal.zt8.2[,4]),log2(chip.signal.zt16.1[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

plot(log2(chip.signal.zt8.2[,4]),log2(chip.signal.zt16.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

chip.signal.zt8.zt16 <- matrix(data = c(chip.signal.zt8.1[,4],
                chip.signal.zt8.2[,4],
                chip.signal.zt16.1[,4],
                chip.signal.zt16.2[,4]),ncol=4)
colnames(chip.signal.zt8.zt16) <- c("ZT8_1","ZT8_2","ZT16_1","ZT16_2")
rownames(chip.signal.zt8.zt16) <- apply(X = chip.signal.zt8.1[,1:3],
                                        MARGIN = 1,
                                        FUN = function(x)
                                          {paste(x,collapse="_")})

boxplot(chip.signal.zt8.zt16,outline=F)

library(FactoMineR)
library(factoextra)
pca.chip.signal.zt8.zt16 <- data.frame(colnames(chip.signal.zt8.zt16),
                                  t(chip.signal.zt8.zt16))
colnames(pca.chip.signal.zt8.zt16)[1] <- "Sample"
res.pca <- PCA(pca.chip.signal.zt8.zt16, 
               graph = FALSE,
               scale.unit = TRUE,
               quali.sup = 1 )

fviz_pca_ind(res.pca, #col.ind = design$group, 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             #addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 2)   

fviz_dend(res.hcpc,k=2,
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 700      # Augment the room for labels
)


chip.signal.zt8.1 <- read.table(
  file = "chip_data/k27_consensus_ll_ld_20/consensus_k27_20_ll_ld/h3k27me3_ld_20_zt8_1_counts_consensus_ll_ld.bed",
  header = F,sep = "\t")
chip.signal.zt8.1[,4] <- chip.signal.zt8.1[,4]/mapped.reads.k27.20c.zt8.1

chip.signal.zt8.2 <- read.table(
  file = "chip_data/k27_consensus_ll_ld_20/consensus_k27_20_ll_ld/h3k27me3_ld_20_zt8_2_counts_consensus_ll_ld.bed",
  header = F,sep = "\t")
chip.signal.zt8.2[,4] <- chip.signal.zt8.2[,4]/mapped.reads.k27.20c.zt8.2

plot(log2(chip.signal.zt8.1[,4]),log2(chip.signal.zt8.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

chip.signal.zt16.1 <- read.table(
  file = "chip_data/k27_consensus_ll_ld_20/consensus_k27_20_ll_ld/h3k27me3_ld_20_zt16_1_counts_consensus_ll_ld.bed",
  header = F,sep = "\t")
chip.signal.zt16.1[,4] <- chip.signal.zt16.1[,4]/mapped.reads.k27.20c.zt16.1

chip.signal.zt16.2 <- read.table(
  file = "chip_data/k27_consensus_ll_ld_20/consensus_k27_20_ll_ld/h3k27me3_ld_20_zt16_2_counts_consensus_ll_ld.bed",
  header = F,sep = "\t")
chip.signal.zt16.2[,4] <- chip.signal.zt16.2[,4]/mapped.reads.k27.20c.zt16.2

plot(log2(chip.signal.zt16.1[,4]),log2(chip.signal.zt16.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

chip.signal.20.1 <- read.table(
  file = "chip_data/k27_consensus_ll_ld_20/consensus_k27_20_ll_ld/h3k27me3_ll_20_1_counts_consensus_ll_ld.bed",
  header = F,sep = "\t")
chip.signal.20.1[,4] <- chip.signal.20.1[,4]/mapped.reads.k27.20c.1

chip.signal.20.2 <- read.table(
  file = "chip_data/k27_consensus_ll_ld_20/consensus_k27_20_ll_ld/h3k27me3_ll_20_2_counts_consensus_ll_ld.bed",
  header = F,sep = "\t")
chip.signal.20.2[,4] <- chip.signal.20.2[,4]/mapped.reads.k27.20c.2

chip.signal.20.3 <- read.table(
  file = "chip_data/k27_consensus_ll_ld_20/consensus_k27_20_ll_ld/h3k27me3_ll_20_3_counts_consensus_ll_ld.bed",
  header = F,sep = "\t")
chip.signal.20.3[,4] <- chip.signal.20.3[,4]/mapped.reads.k27.20c.3


## Differential occupancy analysis

## Changes in H3K27me3 were detected as response to temperature

## All H3K27me3 regions detected at different temperatures are collected in 
## the consensus H3K27me3

h3k27me3.regions <- read.table(file="consensus_ld_ll_all_temps/consensus_k27_ld_20_plus_ll_all_temps/h3k27me3_ld_20_plus_ll_all_temps_consensus.bed",
                    header = F)

h3k27me3.regions <- h3k27me3.regions[,1:3]

nrow(h3k27me3.regions)

#jpeg(filename = "images/h3k27me3_genome_wide_distribution_ll_20.jpg",res = 300,width = 900,height = 900)
#par(mar=c(0,0,0,0))
plot(x=c(0,max.chr.len), y=c(0,4*number.chrs), col = "white", 
     xlab = "", ylab = "", axes=F)  

for(i in 1:number.chrs)
{
  polygon(x = c(0, chr.lens[i], chr.lens[i], 0),   
          y = c(4*i+1, 4*i+1, 4*i-1, 4*i-1),    
          col = "white") 
}

for(i in 1:nrow(h3k27me3.regions))
{
  current.chr <- h3k27me3.regions[i,1]
  current.start <- h3k27me3.regions[i,2]
  current.end <- h3k27me3.regions[i,3]
  
  current.line <- (20:1)[current.chr]
  
  polygon(x = c(current.start, 
                current.end, 
                current.end, 
                current.start),   
          y = c(4*current.line+1, 
                4*current.line+1, 
                4*current.line-1, 
                4*current.line-1),    
          col = "blue",border="blue")
}

h3k27me3.peaks.len <- h3k27me3.regions$V3 - h3k27me3.regions$V2
sum(h3k27me3.peaks.len)/sum(chr.lens)
h3k27me3.regions <- cbind(h3k27me3.regions,h3k27me3.peaks.len)
colnames(h3k27me3.regions) <- c("chr","start","end","length")
head(h3k27me3.regions)
number.of.h3k27.peaks.chr <- vector(mode = "numeric",length = 20)
peaks.len.chr <- vector(mode = "numeric",length = 20)

for(i in 1:20)
{
  number.of.h3k27.peaks.chr[i] <- nrow(subset(h3k27me3.regions, chr == i))
  peaks.len.chr[i] <- sum(subset(h3k27me3.regions, chr == i)[,4]) 
}

percentage.chr.k27 <- peaks.len.chr/rev(chr.lens)

z.scores <- (percentage.chr.k27 - mean(percentage.chr.k27))/sd(percentage.chr.k27)

which(z.scores > 2.5)

p.vals <- pnorm((percentage.chr.k27 - mean(percentage.chr.k27[-19]))/sd(percentage.chr.k27[-19]),lower.tail = F)

p.vals

values.for.colors <- round(-log10(p.vals)*40)
values.for.colors[19] <- 100

white.blue.palette <- colorRampPalette(c('white','blue'))

#jpeg(filename = "images/h3k27me3_chr_barplot_ll_20.jpg",res = 300,width = 1800,height = 900)
barplot(100*percentage.chr.k27,border="black",names.arg = paste0("Chr",1:20),las=2,col=white.blue.palette(100)[values.for.colors]
)
abline(h = 100*mean(percentage.chr.k27[-19]),lty=5)
#dev.off()

library(TxDb.Otauri.JGI)
txdb <- TxDb.Otauri.JGI

ostta.genes <- as.data.frame(genes(txdb))
head(ostta.genes)

target.genes <- vector(mode = "character",length = nrow(h3k27me3.regions))
peak.annotation <- vector(mode = "character",length = nrow(h3k27me3.regions))

for(i in 1:nrow(h3k27me3.regions))
{
  current.chr <- h3k27me3.regions[i,1]
  current.start <- h3k27me3.regions[i,2]
  current.end <- h3k27me3.regions[i,3]
  
  complete.gene.body <- subset(ostta.genes, 
                               seqnames == current.chr & 
                                 start >= current.start & 
                                 end <= current.end)
  
  overlap.tss.pos.strand <- subset(ostta.genes, 
                                   seqnames == current.chr & 
                                     start >= current.start & start <= current.end & 
                                     strand == "+")
  
  overlap.tss.neg.strand <- subset(ostta.genes, 
                                   seqnames == current.chr & 
                                     end >= current.start & end <= current.end & 
                                     strand == "-")
  
  overlap.tes.pos.strand <- subset(ostta.genes, 
                                   seqnames == current.chr & 
                                     end >= current.start & end <= current.end & 
                                     strand == "+")
  
  overlap.tes.neg.strand <- subset(ostta.genes, 
                                   seqnames == current.chr & 
                                     start >= current.start & start <= current.end & 
                                     strand == "-")
  
  inside.gene.body <- subset(ostta.genes, seqnames == current.chr & 
                               start <= current.start & end >= current.end)
  
  
  if(nrow(complete.gene.body) > 0)
  {
    target.genes[i] <- paste(unlist(complete.gene.body$gene_id),
                             collapse = ",")
    peak.annotation[i] <- "gene.body"
  } else if(nrow(overlap.tss.pos.strand) > 0)
  {
    target.genes[i] <- paste(unlist(overlap.tss.pos.strand$gene_id),
                             collapse = ",")
    peak.annotation[i] <- "TSS"
  } else if(nrow(overlap.tss.neg.strand) > 0)
  {
    target.genes[i] <- paste(unlist(overlap.tss.neg.strand$gene_id),
                             collapse = ",")
    peak.annotation[i] <- "TSS"
  } else if(nrow(overlap.tes.pos.strand) > 0)
  {
    target.genes[i] <- paste(unlist(overlap.tes.pos.strand$gene_id),
                             collapse = ",")
    peak.annotation[i] <- "TES"
  } else if(nrow(overlap.tes.neg.strand) > 0)
  {
    target.genes[i] <- paste(unlist(overlap.tes.neg.strand$gene_id),
                             collapse = ",")
    peak.annotation[i] <- "TES"
  } else if(nrow(inside.gene.body) > 0)
  {
    target.genes[i] <- paste(unlist(inside.gene.body$gene_id),
                             collapse = ",")
    peak.annotation[i] <- "internal.gene.body"
  } else
  {
    peak.annotation[i] <- "intergenic"
  }
}

freq.peak.annotation <- table(peak.annotation)
#jpeg(filename = "images/h3k27me3_pie_chart_ll_20.jpg",res = 300,width = 900,height = 900)
#par(mar=c(0,0,0,0))
pie(freq.peak.annotation[c("gene.body",
                           "internal.gene.body",
                           "TSS", 
                           "TES",
                           "intergenic")],
    labels = c("Full Gene Body", 
               "Gene Body", 
               "TSS", 
               "TES", 
               "Intergenic"),
    clockwise = T,border = "black",col=c("blue","cyan","green","red","grey"))
#dev.off()

percentage.freq.peak.annotation <- 
  100*freq.peak.annotation/sum(freq.peak.annotation)


percentage.freq.peak.annotation["gene.body"] +
  percentage.freq.peak.annotation["internal.gene.body"]

percentage.freq.peak.annotation

h3k27me3.regions.peak.annotation <- cbind(h3k27me3.regions[,1:3],
                               data.frame(peak.annotation,target.genes))
nrow(h3k27me3.regions.peak.annotation)

head(h3k27me3.regions.peak.annotation)

write.table(x = h3k27me3.regions.peak.annotation,
            file = "tables/H3K27me3_consensus_peak_annotation.tsv",
            quote = F,sep = "\t",row.names = F)

all.k27.genes <- (unlist(sapply(X =h3k27me3.regions.peak.annotation$target.genes,
                                FUN = function(x) {strsplit(x,",")[[1]]})))
names(all.k27.genes) <- NULL
length(all.k27.genes)
length(all.k27.genes)/7668

dna.transposons <- 
  read.table(file="../te_analysis/dna_transposons_unsorted.bed",
             sep="\t",header=F)
nrow(dna.transposons)

ltr.transposons <- 
  read.table(file="../te_analysis/ltr_transposons_unsorted.bed",
             sep="\t", header=F)
nrow(ltr.transposons)

unclassified.transposons <- 
  read.table(file="../te_analysis/unknown_transposons_unsorted.bed", 
             sep="\t",header = F)
nrow(unclassified.transposons)

overlapping.peaks <- function(target.peak, peaks)
{
  subset(peaks, V1 == target.peak[1] & (
    (V2 >= target.peak[2] & V2 <= target.peak[3]) |
      (V3 >= target.peak[2] & V3 <= target.peak[3]) |
      (V2 <= target.peak[2] & V3 >= target.peak[3])))
}


peak.te <- vector(mode = "character",length = nrow(h3k27me3.regions.peak.annotation))

for(i in 1:nrow(h3k27me3.regions.peak.annotation))
{
  current.peak <- unlist(h3k27me3.regions.peak.annotation[i,1:3])
  
  if(nrow(overlapping.peaks(current.peak,dna.transposons)) > 0)
  {
    peak.te[i] <- "DNA_transposon"
  } else if (nrow(overlapping.peaks(current.peak,ltr.transposons)) > 0)
  {
    peak.te[i] <- "LTR_transposon"
  } else if (nrow(overlapping.peaks(current.peak,unclassified.transposons)) > 0)
  {
    peak.te[i] <- "Unclassified_transposon"
  } else
  {
    peak.te[i] <- "No_transposon"
  }
}

100*table(peak.te)/sum(table(peak.te))

pie(table(peak.te), col=c("black","grey","white","lightblue"))

nrow(dna.transposons)
100*table(peak.te)["DNA_transposon"]/nrow(dna.transposons)

nrow(ltr.transposons)
100*table(peak.te)["LTR_transposon"]/nrow(ltr.transposons)

nrow(unclassified.transposons)
100*table(peak.te)["Unclassified_transposon"]/nrow(unclassified.transposons)


table(h3k27me3.regions.peak.annotation$peak.annotation[peak.te == "DNA_transposon"])
table(h3k27me3.regions.peak.annotation$peak.annotation[peak.te == "LTR_transposon"])
table(h3k27me3.regions.peak.annotation$peak.annotation[peak.te == "Unclassified_transposon"])

table(h3k27me3.regions.peak.annotation$peak.annotation)

(table(h3k27me3.regions.peak.annotation$peak.annotation[peak.te == "DNA_transposon"])["intergenic"] +
    table(h3k27me3.regions.peak.annotation$peak.annotation[peak.te == "LTR_transposon"])["intergenic"] +
    table(h3k27me3.regions.peak.annotation$peak.annotation[peak.te == "Unclassified_transposon"])["intergenic"]) / table(h3k27me3.regions.peak.annotation$peak.annotation)["intergenic"]

h3k27me3.regions.peak.annotation <- cbind(h3k27me3.regions.peak.annotation,peak.te)
write.table(x = h3k27me3.regions.peak.annotation,
            file = "tables/H3K27me3_consensus_peak_annotation.tsv",quote = F,sep = "\t",row.names = F)


library(clusterProfiler)
library(enrichplot)
library(org.Otauriv5.eg.db)
k27.enrichment <- enrichGO(gene = all.k27.genes, pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           OrgDb = org.Otauriv5.eg.db,
                           ont = "BP", 
                           keyType = "GID")
k27.enrichment.df <- as.data.frame(k27.enrichment)

k27.go.terms <- enrichGO(gene = all.k27.genes, pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         OrgDb = org.Otauriv5.eg.db,
                         ont = "BP", 
                         keyType = "GID")
k27.go.terms.filtered <- k27.go.terms
k27.go.terms.filtered@result <- subset(k27.go.terms@result, pvalue < 0.05)

k27.go.terms.df <- as.data.frame(k27.go.terms.filtered)


write.table(x = k27.go.terms.df, file = "go_enrichment_k27.tsv",
            quote = F,sep = "\t",row.names = F)

treeplot(x = pairwise_termsim(k27.go.terms.filtered))


k27.enrichment <- enrichGO(gene = k27.gene.body, pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           OrgDb = org.Otauriv5.eg.db,
                           ont = "BP", 
                           keyType = "GID")
k27.enrichment.df <- as.data.frame(k27.enrichment)

k27.go.terms <- enrichGO(gene = k27.gene.body, pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         OrgDb = org.Otauriv5.eg.db,
                         ont = "BP", 
                         keyType = "GID")
k27.go.terms.filtered <- k27.go.terms
k27.go.terms.filtered@result <- subset(k27.go.terms@result, pvalue < 0.05)

k27.go.terms.df <- as.data.frame(k27.go.terms.filtered)

##jpeg(filename = "../images/h3k27me3_treeplot_gene_body.jpg",res = 250,width = 3800,height = 1500)
treeplot(x = pairwise_termsim(k27.go.terms.filtered))
##dev.off()

gene.ids <- ostta.genes$gene_id
k27.enrich.kegg <- enrichKEGG(gene  = paste0("OT_",all.k27.genes),
                              universe = paste0("OT_",gene.ids),
                              organism = "ota",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05)

df.activated.enrich.kegg <- as.data.frame(k27.enrich.kegg)
head(df.activated.enrich.kegg)

library(pathview)

pathview.genes <- rep(0,length(gene.ids))
names(pathview.genes) <- gene.ids

pathview.genes[all.k27.genes] <- 1

names(pathview.genes) <- paste0("OT_",gene.ids)

pathview(gene.data = sort(pathview.genes,decreasing = TRUE),
         pathway.id = "ota04814",
         species = "ota",
         gene.idtype = "KEGG")

rownames(h3k27me3.regions.peak.annotation) <-
apply(X = h3k27me3.regions.peak.annotation[,1:3],
      MARGIN = 1,
      FUN = function(x){paste(x,collapse = "_")})


mapped.reads.k27.20c.zt8.1 <- 2254676/1E6 
mapped.reads.k27.20c.zt8.2 <- 2639017/1E6
mapped.reads.k27.20c.zt16.1 <- 3605883/1E6
mapped.reads.k27.20c.zt16.2 <- 3405826/1E6

mapped.reads.k27.10c.1 <- 6764067/1E6 
mapped.reads.k27.10c.2 <- 6429853/1E6
mapped.reads.k27.10c.3 <- 6625634/1E6

mapped.reads.k27.14c.1 <- 24537863/1E6 
mapped.reads.k27.14c.2 <- 6592650/1E6
mapped.reads.k27.14c.3 <- 4317544/1E6

mapped.reads.k27.20c.1 <- 3146662/1E6 
mapped.reads.k27.20c.2 <- 7702115/1E6
mapped.reads.k27.20c.3 <- 7102648/1E6

mapped.reads.k27.26c.1 <- 5635098/1E6 
mapped.reads.k27.26c.2 <- 4986577/1E6
mapped.reads.k27.26c.3 <- 3988645/1E6

chip.signal.zt8.1 <- read.table(
  file = "h3k27me3/h3k27me3_ld_20_zt8_chip_1_counts.bed",
  header = F,sep = "\t")
chip.signal.zt8.1[,4] <- chip.signal.zt8.1[,4]/mapped.reads.k27.20c.zt8.1

chip.signal.zt8.2 <- read.table(
  file = "h3k27me3/h3k27me3_ld_20_zt8_chip_2_counts.bed",
  header = F,sep = "\t")
chip.signal.zt8.2[,4] <- chip.signal.zt8.2[,4]/mapped.reads.k27.20c.zt8.2

chip.signal.zt16.1 <- read.table(
  file = "h3k27me3/h3k27me3_ld_20_zt16_chip_1_counts.bed",
  header = F,sep = "\t")
chip.signal.zt16.1[,4] <- chip.signal.zt16.1[,4]/mapped.reads.k27.20c.zt16.1

chip.signal.zt16.2 <- read.table(
  file = "h3k27me3/h3k27me3_ld_20_zt16_chip_2_counts.bed",
  header = F,sep = "\t")
chip.signal.zt16.2[,4] <- chip.signal.zt16.2[,4]/mapped.reads.k27.20c.zt16.2

chip.signal.20.1 <- read.table(
  file = "h3k27me3/h3k27me3_ll_20_chip_1_counts.bed",
  header = F,sep = "\t")
chip.signal.20.1[,4] <- chip.signal.20.1[,4]/mapped.reads.k27.20c.1

chip.signal.20.2 <- read.table(
  file = "h3k27me3/h3k27me3_ll_20_chip_2_counts.bed",
  header = F,sep = "\t")
chip.signal.20.2[,4] <- chip.signal.20.2[,4]/mapped.reads.k27.20c.2

chip.signal.20.3 <- read.table(
  file = "h3k27me3/h3k27me3_ll_20_chip_3_counts.bed",
  header = F,sep = "\t")
chip.signal.20.3[,4] <- chip.signal.20.3[,4]/mapped.reads.k27.20c.3

chip.signal.10.1 <- read.table(
  file = "h3k27me3/h3k27me3_ll_10_chip_1_counts.bed",
  header = F,sep = "\t")
chip.signal.10.1[,4] <- chip.signal.10.1[,4]/mapped.reads.k27.10c.1

chip.signal.10.2 <- read.table(
  file = "h3k27me3/h3k27me3_ll_10_chip_2_counts.bed",
  header = F,sep = "\t")
chip.signal.10.2[,4] <- chip.signal.10.2[,4]/mapped.reads.k27.10c.2

chip.signal.10.3 <- read.table(
  file = "h3k27me3/h3k27me3_ll_10_chip_3_counts.bed",
  header = F,sep = "\t")
chip.signal.10.3[,4] <- chip.signal.10.3[,4]/mapped.reads.k27.10c.3

chip.signal.14.1 <- read.table(
  file = "h3k27me3/h3k27me3_ll_14_chip_1_counts.bed",
  header = F,sep = "\t")
chip.signal.14.1[,4] <- chip.signal.14.1[,4]/mapped.reads.k27.14c.1

chip.signal.14.2 <- read.table(
  file = "h3k27me3/h3k27me3_ll_14_chip_2_counts.bed",
  header = F,sep = "\t")
chip.signal.14.2[,4] <- chip.signal.14.2[,4]/mapped.reads.k27.14c.2

chip.signal.14.3 <- read.table(
  file = "h3k27me3/h3k27me3_ll_14_chip_3_counts.bed",
  header = F,sep = "\t")
chip.signal.14.3[,4] <- chip.signal.14.3[,4]/mapped.reads.k27.14c.3

chip.signal.26.1 <- read.table(
  file = "h3k27me3/h3k27me3_ll_26_chip_1_counts.bed",
  header = F,sep = "\t")
chip.signal.26.1[,4] <- chip.signal.26.1[,4]/mapped.reads.k27.26c.1

chip.signal.26.2 <- read.table(
  file = "h3k27me3/h3k27me3_ll_26_chip_2_counts.bed",
  header = F,sep = "\t")
chip.signal.26.2[,4] <- chip.signal.26.2[,4]/mapped.reads.k27.26c.2

chip.signal.26.3 <- read.table(
  file = "h3k27me3/h3k27me3_ll_26_chip_3_counts.bed",
  header = F,sep = "\t")
chip.signal.26.3[,4] <- chip.signal.26.3[,4]/mapped.reads.k27.26c.3

plot(log2(chip.signal.zt8.1[,4]),log2(chip.signal.zt8.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.zt16.1[,4]),log2(chip.signal.zt16.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")


plot(log2(chip.signal.10.1[,4]),log2(chip.signal.10.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.10.1[,4]),log2(chip.signal.10.3[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.10.2[,4]),log2(chip.signal.10.3[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

plot(log2(chip.signal.14.1[,4]),log2(chip.signal.14.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.14.1[,4]),log2(chip.signal.14.3[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.14.2[,4]),log2(chip.signal.14.3[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

plot(log2(chip.signal.20.1[,4]),log2(chip.signal.20.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.20.1[,4]),log2(chip.signal.20.3[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.20.2[,4]),log2(chip.signal.20.3[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")

plot(log2(chip.signal.26.1[,4]),log2(chip.signal.26.2[,4]))
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.26.1[,4]),log2(chip.signal.26.3[,4])+0.3)
lines(x=c(0,20),y=c(0,20),col="blue")
plot(log2(chip.signal.26.2[,4]),log2(chip.signal.26.3[,4])+log2(2^0.3))
lines(x=c(0,20),y=c(0,20),col="blue")



chip.signal <- matrix(data = c(chip.signal.zt8.1[,4],
                               chip.signal.zt8.2[,4],
                               chip.signal.zt16.1[,4],
                               chip.signal.zt16.2[,4],
                               chip.signal.10.1[,4],
                               chip.signal.10.2[,4],
                               chip.signal.10.3[,4],
                               chip.signal.14.1[,4],
                               chip.signal.14.2[,4],
                               chip.signal.14.3[,4],
                               chip.signal.20.1[,4],
                               chip.signal.20.2[,4],
                               chip.signal.20.3[,4],
                               chip.signal.26.1[,4],
                               chip.signal.26.2[,4],
                               chip.signal.26.3[,4]),ncol=16)
colnames(chip.signal) <- c("ZT8_1","ZT8_2",
                           "ZT16_1","ZT16_2",
                           "t10C_1","t10C_2","t10C_3",
                           "t14C_1","t14C_2","t14C_3",
                           "t20C_1","t20C_2","t20C_3",
                           "t26C_1","t26C_2","t26C_3")
rownames(chip.signal) <- apply(X = chip.signal.10.1[,1:3],
                               MARGIN = 1,
                               FUN = 
                                 function(x)
                                   {
                                     paste(x[1],paste(x[2],x[3],sep="-"),sep=":")
                                   }
                               )

boxplot(log2(chip.signal),outline=F)


head(chip.signal)
chip.signal.zt8 <- rowMeans(chip.signal[,paste("ZT8",1:2,sep="_")])
chip.signal.zt16 <- rowMeans(chip.signal[,paste("ZT16",1:2,sep="_")])
chip.signal.10 <- rowMeans(chip.signal[,paste("t10C",1:3,sep="_")])
chip.signal.14 <- rowMeans(chip.signal[,paste("t14C",1:3,sep="_")])
chip.signal.20 <- rowMeans(chip.signal[,paste("t20C",1:3,sep="_")])
chip.signal.26 <- rowMeans(chip.signal[,paste("t26C",1:3,sep="_")])


plot(log2(chip.signal.zt8),log2(chip.signal.zt16))
lines(x=c(0,20),y=c(0,20),lwd=3,col="blue")

plot(log2(chip.signal.10),log2(chip.signal.14))
lines(x=c(0,20),y=c(0,20),lwd=3,col="blue")

plot(log2(chip.signal.10),log2(chip.signal.20))
lines(x=c(0,20),y=c(0,20),lwd=3,col="blue")

plot(log2(chip.signal.10),log2(chip.signal.26))
lines(x=c(0,20),y=c(0,20),lwd=3,col="blue")


library(FactoMineR)
library(factoextra)

pca.temp.chip.signal <- data.frame(colnames(normalized.temp.chip.signal),
                                   t(normalized.temp.chip.signal))


pca.temp.chip.signal <- data.frame(colnames(chip.signal),
                                   t(chip.signal))

colnames(pca.temp.chip.signal)[1] <- "Sample"
res.pca <- PCA(pca.temp.chip.signal, 
               graph = FALSE,
               scale.unit = TRUE,
               quali.sup = 1 )

fviz_pca_ind(res.pca, col.ind = c("LD","LD","LD","LD",
                                  rep("LL10",3),rep("LL14",3),
                                  rep("LL20",3),rep("LL26",3)),#exp.design$group, 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 4)   

fviz_dend(res.hcpc,k=4,
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 50      # Augment the room for labels
)



log2.chip.signal <- log2(chip.signal)

log2.chip.signal <- normalized.temp.chip.signal
log2.chip.signal["8:624229-627058",]
chip.signal["8:624229-627058",]


t10.t26["8:624229-627058",]
log2(1.25)


library(limma)
limma.experimental.design <- model.matrix(~ -1+#factor(c(1,1,1,
                                            #            2,2,2,
                                            #            3,3,3,
                                            #            4,4,4)))
                                            factor(c(1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6)))
colnames(limma.experimental.design) <- c("LDZT8","LDZT16",
                                         "LL10","LL14","LL20","LL26")

linear.fit <- lmFit(log2.chip.signal, limma.experimental.design)

contrast.matrix <- makeContrasts(LDZT16-LDZT8,
                                 LL26-LL10,
                                 LL20-LL10,
                                 LL14-LL10, 
                                 levels = c("LDZT8","LDZT16","LL10","LL14","LL20","LL26"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)                               

t10.t26 <- topTable(contrast.results, 
                     number=nrow(log2.chip.signal), 
                     coef = 2, sort.by = "logFC", )
head(t10.t26)
t10.t26["8:624229-627058",] #ostta08g03710
t10.t26["18:189253-192477",] #ostta18g01040
t10.t26["16:170969-173026", ] #ostta16g01100
t10.t26["17:52909-54964", ] #ostta17g00300



t10.t20.logfc <- t10.t20$logFC
t10.t20.adj.p.val <- t10.t20$adj.P.Val

higher.peaks.t10.t20 <- rownames(log2.chip.signal)[t10.t20.logfc > log2(1.5) & 
                                                     t10.t20.adj.p.val < 0.05]
length(higher.peaks.t10.t20)

lower.peaks.t10.t20 <- rownames(log2.chip.signal)[t10.t20.logfc < -log2(1.5) & 
                                                    t10.t20.adj.p.val < 0.05]
length(lower.peaks.t10.t20)


t10.t14 <- topTable(contrast.results, 
                     number=nrow(log2.chip.signal), 
                     coef = 2, sort.by = "logFC", )
head(t10.t14)
t10.t14.logfc <- t10.t14$logFC
t10.t14.adj.p.val <- t10.t14$adj.P.Val
plot(t10.t14.logfc,-log10(t10.t14.adj.p.val))

higher.peaks.t10.t14.1.25 <- rownames(log2.chip.signal)[t10.t14.logfc > log2(1.25) & 
                                                 t10.t14.adj.p.val < 0.05]
higher.peaks.t10.t14.1.5 <- rownames(log2.chip.signal)[t10.t14.logfc > log2(1.5) & 
                                                          t10.t14.adj.p.val < 0.05]
higher.peaks.t10.t14.1.75 <- rownames(log2.chip.signal)[t10.t14.logfc > log2(1.75) & 
                                                         t10.t14.adj.p.val < 0.05]
higher.peaks.t10.t14.2 <- rownames(log2.chip.signal)[t10.t14.logfc > log2(2) & 
                                                          t10.t14.adj.p.val < 0.05]

length(higher.peaks.t10.t14)
lower.t10.t14 <- rownames(log2.chip.signal)[t10.t14.logfc < -1 & 
                                                 t10.t14.adj.p.val < 0.05]
length(lower.t10.t14)



t10.t20 <- topTable(contrast.results, 
                    number=nrow(log2.chip.signal), 
                    coef = 3, sort.by = "logFC", )
head(t10.t20)
t10.t20.logfc <- t10.t20$logFC
t10.t20.adj.p.val <- t10.t20$adj.P.Val

plot(t10.t20.logfc,-log10(t10.t20.adj.p.val))

higher.peaks.t10.t20.1.25 <- rownames(log2.chip.signal)[t10.t20.logfc > log2(1.25) & 
                                                     t10.t20.adj.p.val < 0.05]
higher.peaks.t10.t20.1.5 <- rownames(log2.chip.signal)[t10.t20.logfc > log2(1.5) & 
                                                          t10.t20.adj.p.val < 0.05]
higher.peaks.t10.t20.1.75 <- rownames(log2.chip.signal)[t10.t20.logfc > log2(1.75) & 
                                                         t10.t20.adj.p.val < 0.05]
higher.peaks.t10.t20.2 <- rownames(log2.chip.signal)[t10.t20.logfc > log2(2) & 
                                                          t10.t20.adj.p.val < 0.05]

length(higher.peaks.t10.t20)
lower.t10.t20 <- rownames(log2.chip.signal)[t10.t20.logfc < -1 & 
                                              t10.t20.adj.p.val < 0.05]
length(lower.t10.t20)



t10.t26 <- topTable(contrast.results, 
                    number=nrow(log2.chip.signal), 
                    coef = 4, sort.by = "logFC", )
head(t10.t26)
t10.t26["18:189254-192477",]
t10.t26.logfc <- t10.t26$logFC
t10.t26.adj.p.val <- t10.t26$adj.P.Val
plot(t10.t26.logfc, -log10(t10.t26.adj.p.val))
higher.peaks.t10.t26.1.25 <- rownames(log2.chip.signal)[t10.t26.logfc > log2(1.25) & 
                                                     t10.t26.adj.p.val < 0.05]
higher.peaks.t10.t26.1.5 <- rownames(log2.chip.signal)[t10.t26.logfc > log2(1.5) & 
                                                          t10.t26.adj.p.val < 0.05]
higher.peaks.t10.t26.1.75 <- rownames(log2.chip.signal)[t10.t26.logfc > log2(1.75) & 
                                                          t10.t26.adj.p.val < 0.05]
higher.peaks.t10.t26.2 <- rownames(log2.chip.signal)[t10.t26.logfc > log2(2) & 
                                                          t10.t26.adj.p.val < 0.05]

length(higher.peaks.t10.t26.1.25)
lower.t10.t26 <- rownames(log2.chip.signal)[t10.t26.logfc < -1 & 
                                              t10.t26.adj.p.val < 0.05]
length(lower.t10.t26)





boxplot(chip.signal[higher.peaks.t10.t26.2,],outline=F)


boxplot(chip.signal.10[higher.peaks.t10.t26.1.25],
chip.signal.14[higher.peaks.t10.t26.1.25],
chip.signal.20[higher.peaks.t10.t26.1.25],
chip.signal.26[higher.peaks.t10.t26.1.25],outline=F)

peak.igv <- function(x)
{
  paste(x[1],paste(x[2],x[3],sep="-"),sep=":")
}

chip.signal[apply(X = subset(ll.20.peak.annotation, target.genes %in% activated.genes.10.vs.20)[,1:3],MARGIN = 1,FUN = peak.igv),]

summary(t10.t26[apply(X = subset(ll.20.peak.annotation, target.genes %in% activated.genes.10.vs.20)[,1:3],MARGIN = 1,FUN = peak.igv),"logFC"])
2^0.6
boxplot(t10.t26[apply(X = subset(ll.20.peak.annotation, target.genes %in% activated.genes.10.vs.20)[,1:3],MARGIN = 1,FUN = peak.igv),"logFC"],
        t10.t26$logFC)


genes.increased.t10.t20 <- h3k27me3.regions.peak.annotation[activated.t10.t20,"target.genes"]
genes.increased.t10.t20 <- unlist(strsplit(x = genes.increased.t10.t20,split=","))

head(gene.expression)

normalized.gene.expression <- 2^read.table(file = "/home/fran/Dropbox/aveiro/group0/ostta_temp/Quantile-normalized.txt",header = T)-1
head(normalized.gene.expression)
boxplot(normalized.gene.expression,outline=F)

gene.ids <- read.table(file = "/home/fran/Dropbox/aveiro/group0/ostta_26C_100uE_3.tsv",header=T,sep="\t")$Gene.ID

rownames(normalized.gene.expression) <- gene.ids

write.table(x = normalized.gene.expression,file = "normalized_gene_expression.tsv",quote = F,sep = "\t")

normalized.gene.expression <- read.table(file = "normalized_gene_expression.tsv",header = T,sep = "\t")
head(normalized.gene.expression)
normalized.gene.expression <- gene.expression


gene.id <- "ostta17g00300"
gene.id <- "ostta02g00720"
gene.id <- "ostta17g00300" #ok
gene.id <- "ostta03g02570"
gene.id <- "ostta02g00780"
gene.id <- "ostta18g01040"
gene.id <- "ostta08g03710" #ok
gene.id <- "ostta06g04460"
gene.id <- "ostta14g00070"
gene.id <- "ostta14g02090" #ok AP2
gene.id <- "ostta17g00450" #ok BELL1
gene.id <- "ostta11g02490" #ok WRKY
gene.id <- "ostta07g04340"
gene.id <- "ostta02g03030" #ok CSP
gene.id <- "ostta15g01070"
gene.id <- "ostta09g00590"
gene.id <- "ostta02g04200"
gene.expression.10 <- unlist(normalized.gene.expression[gene.id,paste("t10C",1:3,sep="_")])
gene.expression.14 <- unlist(normalized.gene.expression[gene.id,paste("t14C",1:3,sep="_")])
gene.expression.20 <- unlist(normalized.gene.expression[gene.id,paste("t20C",1:3,sep="_")])
gene.expression.26 <- unlist(normalized.gene.expression[gene.id,paste("t26C",1:3,sep="_")])

means <- c(mean(gene.expression.10),
           mean(gene.expression.14),
           mean(gene.expression.20),
           mean(gene.expression.26))

sds <- c(sd(gene.expression.10),
         sd(gene.expression.14),
         sd(gene.expression.20),
         sd(gene.expression.26))

y.max <- 1.2*max(means + sds)

par(lwd=2)
xpos <- barplot(means,ylim=c(0,y.max),
        col = colorRampPalette(c("blue", "red"))(4),cex.axis = 1.5,
        names.arg = c("10ºC","14ºC","20ºC","26ºC"),cex.names=2,lwd=2)
par(lwd=1)
arrows(x0 = xpos,y0 = means-sds/sqrt(3), x1 = xpos,y1=means+sds/sqrt(3),code=3,angle=90,length=0.07,lwd=2,border="black")


points(x = rep(xpos+0.2,each=3),y=c(gene.expression.10,
                                gene.expression.14,
                                gene.expression.20,
                                gene.expression.26),lwd=2)









activated.enrich.go <- enrichGO(gene          = intersect(all.k27.genes,activated.genes.10.vs.20),
                                OrgDb         = org.Otauriv5.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                readable      = FALSE,
                                keyType = "GID")

as.data.frame(activated.enrich.go)

rownames(ll.20.peak.annotation) <- rownames(chip.signal)

higher.peaks.genes <- ll.20.peak.annotation[higher.peaks.t10.t26.2,"target.genes"]


boxplot(gene.expression.10[higher.peaks.genes],
        gene.expression.14[higher.peaks.genes],
        gene.expression.20[higher.peaks.genes],
        gene.expression.26[higher.peaks.genes],outline=F)


genes.increased.t10.t20 <- higher.peaks.genes

library(clusterProfiler)
library(org.Otauriv5.eg.db)
k27.enrichment <- enrichGO(gene = genes.increased.t10.t20, pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           OrgDb = org.Otauriv5.eg.db,
                           ont = "BP", 
                           keyType = "GID")
k27.enrichment.df <- as.data.frame(k27.enrichment)

k27.go.terms <- enrichGO(gene = genes.increased.t10.t20, pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         OrgDb = org.Otauriv5.eg.db,
                         ont = "BP", 
                         keyType = "GID")
k27.go.terms.filtered <- k27.go.terms
k27.go.terms.filtered@result <- subset(k27.go.terms@result, pvalue < 0.05)

k27.go.terms.df <- as.data.frame(k27.go.terms.filtered)

##jpeg(filename = "../images/h3k27me3_treeplot_gene_body.jpg",res = 250,width = 3800,height = 1500)
treeplot(x = pairwise_termsim(k27.go.terms.filtered))

gene.ids <- ostta.genes$gene_id
k27.enrich.kegg <- enrichKEGG(gene  = paste0("OT_",genes.increased.t10.t20),
                              universe = paste0("OT_",gene.ids),
                              organism = "ota",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05)

df.activated.enrich.kegg <- as.data.frame(k27.enrich.kegg)
head(df.activated.enrich.kegg)

library(pathview)

pathview.genes <- rep(0,length(gene.ids))
names(pathview.genes) <- gene.ids

pathview.genes[activated.genes.10.vs.20] <- 1

names(pathview.genes) <- paste0("OT_",gene.ids)

pathview(gene.data = sort(pathview.genes,decreasing = TRUE),
         pathway.id = "ota04814",
         species = "ota",
         gene.idtype = "KEGG")


gene.expression.barplot(gene="ostta14g00070", expression.matrix) #myosin
gene.expression.barplot(gene="ostta18g01040", expression.matrix) #myosin
gene.expression.barplot(gene="ostta06g03320", expression.matrix) #myosin
gene.expression.barplot(gene="ostta06g04460", expression.matrix) #myosin
gene.expression.barplot(gene="ostta08g03710", expression.matrix) #myosin

k27.enrichment <- enrichGO(gene = genes.increased.t10.t20, pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           OrgDb = org.Otauriv5.eg.db,
                           ont = "BP", 
                           keyType = "GID")
k27.enrichment.df <- as.data.frame(k27.enrichment)





head(chip.signal)
boxplot(chip.signal,outline=F,las=2)
temp.chip.signal <- chip.signal[,c("t10C_1","t10C_2","t10C_3",
                                   "t14C_1","t14C_2","t14C_3",
                                   "t20C_1","t20C_2","t20C_3",
                                   "t26C_1","t26C_2","t26C_3")]
boxplot(temp.chip.signal,outline=F,las=2)

write.table(x = chip.signal+1,file = "chip_signal.tsv",
            quote = F,row.names = F,
            sep = "\t")


exp.design <- data.frame(sample=colnames(chip.signal),
                         group=c(rep("zt8",2),
                                 rep("zt16",2),
                                 rep("t10C",3),
                                 rep("t14C",3),
                                 rep("t20C",3),
                                 rep("t26C",3)))

write.table(x = exp.design,file = "temp_normalyzer_design.tsv",quote = F,row.names = F,
            sep = "\t")


library(NormalyzerDE)

normalyzer(jobName = "norm_chip",designPath = "temp_normalyzer_design.tsv",
           dataPath = "chip_signal.tsv",outputDir = ".")


normalized.temp.chip.signal <- read.table(file="norm_chip/Quantile-normalized.txt", 
                                         header=T)
head(normalized.temp.chip.signal)
rownames(normalized.temp.chip.signal) <- rownames(temp.chip.signal)
boxplot(normalized.temp.chip.signal)


library(FactoMineR)
library(factoextra)

pca.temp.chip.signal <- data.frame(colnames(normalized.temp.chip.signal),
                                       t(normalized.temp.chip.signal))


pca.temp.chip.signal <- data.frame(colnames(chip.signal),
                                   t(chip.signal))

colnames(pca.temp.chip.signal)[1] <- "Sample"
res.pca <- PCA(pca.temp.chip.signal, 
               graph = FALSE,
               scale.unit = TRUE,
               quali.sup = 1 )

fviz_pca_ind(res.pca, col.ind = c("LD","LD","LD","LD",
                                  rep("LL10",3),rep("LL14",3),
                                  rep("LL20",3),rep("LL26",3)),#exp.design$group, 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 4)   

fviz_dend(res.hcpc,k=4,
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 50      # Augment the room for labels
)











ostta.tfs <- read.table(file = "transcription_factors_list.tsv",header = T,as.is=T)

intersect(ostta.tfs$gene_id,all.k27.genes)
intersect(ostta.tfs$gene_id,k27.gene.body)
intersect(ostta.tfs$gene_id,k27.tss)
intersect(ostta.tfs$gene_id,k27.tes)

subset(ostta.tfs, gene_id %in% all.k27.genes)


genes.data.df <- as.data.frame(genes(txdb))
cds.data <- as.data.frame(cds(txdb))
exons.data <- as.data.frame(exons(txdb))
promoter_length <- 1000
gene.name <- "ostta02g03030"
selected.bigwig.files <- "h3k27me3/h3k27me3_ll_20_chip_1.bw"
selected.bed.files <- "h3k27me3/h3k27me3_ll_20_peaks.bed"

gene.profile <- function(gene.name, 
                         promoter_length, 
                         selected.bigwig.files, 
                         selected.bed.files)
{
  target.gene.body <- genes.data.df[gene.name,]
  target.gene.chr <- as.character(target.gene.body$seqnames)
  target.gene.start <- target.gene.body$start
  target.gene.end <- target.gene.body$end
  
  target.gene.strand <- as.character(target.gene.body$strand)
  
  ## Extract cds annotation
  cds.data.target.gene <- subset(cds.data, 
                                 seqnames == target.gene.chr & 
                                   (start >= target.gene.start & 
                                      end <= target.gene.end))
  
  ## Extract exons annotation
  exons.data.target.gene <- subset(exons.data, 
                                   seqnames == target.gene.chr & 
                                     (start >= target.gene.start & 
                                        end <= target.gene.end))
  
  ## Determine the genome range to plot including promoter, gene body and 5' UTR
  ## This depends on whether the gene is on the forward or reverse strand
  range.to.plot <- target.gene.body
  
  if(target.gene.strand == "+")
  {
    range.to.plot$start <- range.to.plot$start - promoter_length
    range.to.plot$end <- range.to.plot$end + promoter_length
  } else if (target.gene.strand == "-")
  {
    range.to.plot$end <- range.to.plot$end + promoter_length
    range.to.plot$start <- range.to.plot$start - promoter_length
  }
  
  ## Compute the length of the genome range to represent
  current.length <- range.to.plot$end - range.to.plot$start
  
  ## Compute profile in gene
  library(rtracklayer)  
  library(ChIPpeakAnno)
  ## Since ChIPpeakAnno needs more than one region to plot our region
  ## is duplicated 
  regions.plot <- GRanges(rbind(range.to.plot,range.to.plot))
  
  ## Import signal from the bigwig files
  cvglists <- sapply(selected.bigwig.files, import, 
                     format="BigWig", 
                     which=regions.plot, 
                     as="RleList")
  
  ## Compute signal in the region to plot
  chip.signal <- featureAlignedSignal(cvglists, regions.plot, 
                                      upstream=ceiling(current.length/2), 
                                      downstream=ceiling(current.length/2),
                                      n.tile=current.length) 
  
  ## Compute mean signal 
  if(target.gene.strand == "+")
  {
    chip.signal.mean <- colMeans(chip.signal[[1]],na.rm = TRUE)
  } else if (target.gene.strand == "-")
  {
    chip.signal.mean <- rev(colMeans(chip.signal[[1]],na.rm = TRUE))
  }
  
  ## Normalization
  #chip.signal.mean <- 20 * chip.signal.mean / max(chip.signal.mean)
  upper.lim <- 1.1*max(chip.signal.mean)
  ## Colors to draw signal
  line.colors <- "blue"
  area.colors <- "lightblue"
  
  
  cord.x <- 1:current.length
  
  ## Exon width to plot
  exon.width <- 2#1#2
  
  ## Cds width to plot
  cds.width <- 3#1.5#3
  
  ## Width of the rectangule representing the peak reagion
  peak.width <- 1#0.5#1
  
  ## Extract exons for target gene
  exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
  
  ## Transform exon coordinates to current range
  min.pos <- min(exons.data.target.gene$start)
  
  exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + promoter_length
  exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + promoter_length
  
  ## Extract cds for target gene
  cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
  
  cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + promoter_length
  cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + promoter_length
  
  plot(cord.x, rep(gene.height,length(cord.x)),type="l",col="black",lwd=3,ylab="",
       cex.lab=2,axes=FALSE,xlab="",main="",cex.main=2,
       ylim=c(-35,upper.lim),xlim=c(0,max(cord.x)))
  
  ## Represent exons
  for(i in 1:nrow(exons.data.target.gene))
  {
    # Determine start/end for each exon
    current.exon.start <- exons.data.target.gene$start[i]
    current.exon.end <- exons.data.target.gene$end[i]
    
    ## Determine coordinates for each exon polygon and represent it
    exon.x <- c(current.exon.start,current.exon.end,current.exon.end,current.exon.start)
    exon.y <- c(gene.height + exon.width, gene.height + exon.width, gene.height - exon.width, gene.height - exon.width)
    
    polygon(x = exon.x, y = exon.y, col = "blue",border = "blue")
  }
  
  if(nrow(cds.data.target.gene) > 0)
  {
    for(i in 1:nrow(cds.data.target.gene))
    {
      # Determine current cds start/end
      current.cds.start <- cds.data.target.gene$start[i]
      current.cds.end <- cds.data.target.gene$end[i]
      
      # Determine curret cds coordinates for the polygon and represent it
      cds.x <- c(current.cds.start,current.cds.end,current.cds.end,current.cds.start)
      cds.y <- c(gene.height + cds.width, gene.height + cds.width, gene.height - cds.width, gene.height - cds.width)
      polygon(x = cds.x, y = cds.y, col = "blue",border = "blue")
    }
  }
  
  ## Draw arrow to represent transcription direction 
  if(target.gene.strand == "+")
  {
    lines(c(promoter_length,promoter_length,promoter_length+100),
          y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
    lines(c(promoter_length+50,promoter_length+100),y=c(gene.height+6,gene.height+5),lwd=3)
    lines(c(promoter_length+50,promoter_length+100),y=c(gene.height+4,gene.height+5),lwd=3)
  } else if (target.gene.strand == "-")
  {
    lines(c(current.length - promoter_length, current.length - promoter_length, current.length - promoter_length-100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
    lines(c(current.length - promoter_length-50, current.length - promoter_length - 100),y=c(
    gene.height + 6, gene.height + 5),lwd=3)
    lines(c(current.length - promoter_length-50, current.length - promoter_length - 100),y=c(gene.height + 4, gene.height + 5),lwd=3)
  }
  
  ## Draw promoter range
  if(target.gene.strand == "+")
  {
    axis(side = 1,labels = c(- promoter_length, - promoter_length / 2,"TSS"),at = c(1,promoter_length/2,promoter_length),lwd=2,cex=1.5,las=2,cex=2)
  } else if(target.gene.strand == "-")
  {
    axis(side = 1,labels = c("TSS",- promoter_length / 2,- promoter_length),at = c(current.length-promoter_length,current.length-promoter_length/2, current.length),lwd=2,cex=1.5,las=2,cex=2)
  }
  
  ## Draw gene name
  text(x = current.length / 2, y = -33 , 
       labels = bquote(italic(.(gene.name))),cex = 1.7,font = 3)
  
  ## Extract bed file name 1 and read it
  current.peaks <- read.table(file=selected.bed.files,header = F, as.is = T)
  current.peaks <- current.peaks[,1:3]
  colnames(current.peaks) <- c("seqnames","start","end")
  peak.coordinates <- subset(current.peaks, seqnames == as.character(range.to.plot$seqnames) 
                             & start >= range.to.plot$start & end <= range.to.plot$end) 
  current.peaks.to.plot <- peak.coordinates[,2:3]
  
  ## Transform coordinates 
  current.peaks.to.plot <- current.peaks.to.plot - range.to.plot$start
  
  ## Check if there are peaks for the target gene
  if(nrow(current.peaks.to.plot) > 0)
  {
    #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
    for(j in 1:nrow(current.peaks.to.plot))
    {
      ## Extract start and end point of each peak region
      current.peak.start <- current.peaks.to.plot[j,1]
      current.peak.end <- current.peaks.to.plot[j,2]
      
      ## Computer coordinates for polygon and draw it
      peak.x <- c(current.peak.start,current.peak.end,
                  current.peak.end,current.peak.start)
      peak.y <- c(peak.width - 12,   peak.width - 12, 
                  - peak.width - 12, - peak.width - 12)  
      
      polygon(x = peak.x, y = peak.y, col = "blue", border = "blue",lwd=2)
    }
  }
  
  ## Draw profiles if bw file is provided
  ## Compute base line for current TF
  current.base.line <- - 10
  
  ## Represent signal from the current TF
  lines(chip.signal.mean+current.base.line,type="l",col=line.colors,lwd=3)
  ## Determine polygon coordinates and represent it
  cord.y <- c(current.base.line,chip.signal.mean+current.base.line,current.base.line)
  cord.x <- 1:length(cord.y)
  polygon(cord.x,cord.y,col=area.colors)
  
  axis(side = 2,lwd = 2,
       at = seq(from=0,to=ceiling(upper.lim/10)*10,by=10)+current.base.line,
       labels = seq(from=0,to=ceiling(upper.lim/10)*10,by=10))
}


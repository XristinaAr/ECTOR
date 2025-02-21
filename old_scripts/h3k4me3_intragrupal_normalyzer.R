
#Analysis intragrupal with NormalyzerDE (same as H3K4me3.v3)
#Check the same for H3K27me3


mapped.reads.k4.10c.1 <- 4660652/1E6 
mapped.reads.k4.10c.2 <- 4572532/1E6
mapped.reads.k4.10c.3 <- 4490999/1E6

mapped.reads.k4.14c.1 <- 4609469/1E6 
mapped.reads.k4.14c.2 <- 5034594/1E6
mapped.reads.k4.14c.3 <- 5174666/1E6

mapped.reads.k4.20c.1 <- 6832850/1E6 
mapped.reads.k4.20c.2 <- 5993486/1E6
mapped.reads.k4.20c.3 <- 6701422/1E6

mapped.reads.k4.26c.1 <- 6059862/1E6 
mapped.reads.k4.26c.2 <- 5263311/1E6
mapped.reads.k4.26c.3 <- 4251439/1E6

read.counts.k4.temp10c.1 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp10c/h3k4me3_10_chip_1_counts.bed", header = F, sep = "\t" )
read.counts.k4.temp10c.2 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp10c/h3k4me3_10_chip_2_counts.bed", header = F, sep = "\t" )
read.counts.k4.temp10c.3 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp10c/h3k4me3_10_chip_3_counts.bed", header = F, sep = "\t" )

read.counts.k4.temp14c.1 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp14c/h3k4me3_14_chip_1_counts.bed", header = F, sep = "\t" )
read.counts.k4.temp14c.2 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp14c/h3k4me3_14_chip_2_counts.bed", header = F, sep = "\t" )
read.counts.k4.temp14c.3 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp14c/h3k4me3_14_chip_3_counts.bed", header = F, sep = "\t" )

read.counts.k4.temp20c.1 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp20c/h3k4me3_20_chip_1_counts.bed", header = F, sep = "\t" )
read.counts.k4.temp20c.2 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp20c/h3k4me3_20_chip_2_counts.bed", header = F, sep = "\t" )
read.counts.k4.temp20c.3 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp20c/h3k4me3_20_chip_3_counts.bed", header = F, sep = "\t" )


read.counts.k4.temp26c.1 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp26c/h3k4me3_26_chip_1_counts.bed", header = F, sep = "\t" )
read.counts.k4.temp26c.2 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp26c/h3k4me3_26_chip_2_counts.bed", header = F, sep = "\t" )
read.counts.k4.temp26c.3 <- read.table(file = "chipseq_data/bed_data/h3k4me3/temp26c/h3k4me3_26_chip_3_counts.bed", header = F, sep = "\t" )



peak.cpm <- matrix(data = NA, nrow = nrow(read.counts.k4.temp10c.1), ncol = 15 )
colnames(peak.cpm) <-c("chr","start","end","k4_temp10_1","k4_temp10_2","k4_temp10_3", 
                       "k4_temp14_1","k4_temp14_2","k4_temp14_3",
                       "k4_temp20_1","k4_temp20_2","k4_temp20_3",
                       "k4_temp26_1","k4_temp26_2","k4_temp26_3")

peak.cpm[,1]<- read.counts.k4.temp10c.1[,1]
peak.cpm[,2]<- read.counts.k4.temp10c.1[,2]
peak.cpm[,3]<- read.counts.k4.temp10c.1[,3]

peak.cpm[,4]<- read.counts.k4.temp10c.1[,4]/mapped.reads.k4.10c.1
peak.cpm[,5]<- read.counts.k4.temp10c.2[,4]/mapped.reads.k4.10c.2
peak.cpm[,6]<- read.counts.k4.temp10c.3[,4]/mapped.reads.k4.10c.3

peak.cpm[,7]<- read.counts.k4.temp14c.1[,4]/mapped.reads.k4.14c.1
peak.cpm[,8]<- read.counts.k4.temp14c.2[,4]/mapped.reads.k4.14c.2
peak.cpm[,9]<- read.counts.k4.temp14c.3[,4]/mapped.reads.k4.14c.3

peak.cpm[,10]<- read.counts.k4.temp20c.1[,4]/mapped.reads.k4.20c.1
peak.cpm[,11]<- read.counts.k4.temp20c.2[,4]/mapped.reads.k4.20c.2
peak.cpm[,12]<- read.counts.k4.temp20c.3[,4]/mapped.reads.k4.20c.3

peak.cpm[,13]<- read.counts.k4.temp26c.1[,4]/mapped.reads.k4.26c.1
peak.cpm[,14]<- read.counts.k4.temp26c.2[,4]/mapped.reads.k4.26c.2
peak.cpm[,15]<- read.counts.k4.temp26c.3[,4]/mapped.reads.k4.26c.3

head(peak.cpm)

peak.cpm <- as.data.frame(peak.cpm)



boxplot((peak.cpm$k4_temp10_1 + peak.cpm$k4_temp10_2 +peak.cpm$k4_temp10_3)/3,
        (peak.cpm$k4_temp14_1 + peak.cpm$k4_temp14_2 +peak.cpm$k4_temp14_3)/3,
        (peak.cpm$k4_temp20_1 + peak.cpm$k4_temp20_2 +peak.cpm$k4_temp20_3)/3,
        (peak.cpm$k4_temp26_1 + peak.cpm$k4_temp26_2+peak.cpm$k4_temp26_3)/3,
        outline=F,col=c("#00FFFF","#6699FF","#CC33FF","#FF0000"),ylab="CPM",names=c("10ºC","14ºC","20ºC","26ºC"),las=2,cex.lab=2)




library(NormalyzerDE)

#temp10

write.table(x = peak.cpm[,4:6],file = "peak_cpm_temp10.tsv" ,
            quote = F,row.names = F,
            sep = "\t")

design.10 <- data.frame(sample= paste("k4_temp10", 1:3, sep = "_"),
                        group=c(rep("temp10",3)))

write.table(x = design.10,file = "normalyzer_design_k4_temp10.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "temp10_k4_normalization",designPath = "normalyzer_design_k4_temp10.tsv",
           dataPath = "peak_cpm_temp10.tsv",outputDir = ".", requireReplicates = FALSE)


normalized.peak.cpm.temp10 <- read.table(file="temp10_k4_normalization/Quantile-normalized.txt", header=T)
head(normalized.peak.cpm.temp10)


boxplot(normalized.peak.cpm.temp10, outline=F,col=rainbow(3),
        cex.lab=1.5)

#temp14

write.table(x = peak.cpm[,7:9],file = "peak_cpm_temp14.tsv" ,
            quote = F,row.names = F,
            sep = "\t")

design.14 <- data.frame(sample= paste("k4_temp14", 1:3, sep = "_"),
                        group=c(rep("temp14",3)))

write.table(x = design.14,file = "normalyzer_design_k4_temp14.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "temp14_k4_normalization",designPath = "normalyzer_design_k4_temp14.tsv",
           dataPath = "peak_cpm_temp14.tsv",outputDir = ".", requireReplicates = FALSE)


normalized.peak.cpm.temp14 <- read.table(file="temp14_k4_normalization/Quantile-normalized.txt", header=T)
head(normalized.peak.cpm.temp14)


boxplot(normalized.peak.cpm.temp14, outline=F,col=rainbow(3),
        cex.lab=1.5)

#temp20

write.table(x = peak.cpm[,10:12],file = "peak_cpm_temp20.tsv" ,
            quote = F,row.names = F,
            sep = "\t")

design.20 <- data.frame(sample= paste("k4_temp20", 1:3, sep = "_"),
                        group=c(rep("temp20",3)))

write.table(x = design.20,file = "normalyzer_design_k4_temp20.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "temp20_k4_normalization",designPath = "normalyzer_design_k4_temp20.tsv",
           dataPath = "peak_cpm_temp20.tsv",outputDir = ".", requireReplicates = FALSE)


normalized.peak.cpm.temp20 <- read.table(file="temp20_k4_normalization/Quantile-normalized.txt", header=T)
head(normalized.peak.cpm.temp20)


boxplot(normalized.peak.cpm.temp20, outline=F,col=rainbow(3),
        cex.lab=1.5)


#temp26

write.table(x = peak.cpm[,16:18],file = "peak_cpm_temp26.tsv" ,
            quote = F,row.names = F,
            sep = "\t")


design.26 <- data.frame(sample= paste("k4_temp26", 1:3, sep = "_"),
                        group=c(rep("temp26",3)))

write.table(x = design.26,file = "normalyzer_design_k4_temp26.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "temp26_k4_normalization",designPath = "normalyzer_design_k4_temp26.tsv",
           dataPath = "peak_cpm_temp26.tsv",outputDir = ".", requireReplicates = FALSE)


normalized.peak.cpm.temp26 <- read.table(file="temp26_k4_normalization/Quantile-normalized.txt", header=T)
head(normalized.peak.cpm.temp26)



boxplot(normalized.peak.cpm.temp26, outline=F,col=rainbow(4),
        cex.lab=1.5)





normalized.peak.cpm <- cbind(normalized.peak.cpm.temp10,normalized.peak.cpm.temp14,
                             normalized.peak.cpm.temp20,
                             normalized.peak.cpm.temp26)
head(normalized.peak.cpm)

peak.cpm[,1]<- read.counts.k4.temp10c.1[,1]
peak.cpm[,2]<- read.counts.k4.temp10c.1[,2]
peak.cpm[,3]<- read.counts.k4.temp10c.1[,3]


normalized.peak.cpm <-  cbind(peak.cpm[,1], peak.cpm[,2], peak.cpm[,3], normalized.peak.cpm)
colnames(normalized.peak.cpm)[1:3] <- c("chr", "start", "end")
head(normalized.peak.cpm)



boxplot((normalized.peak.cpm$k4_temp10_1 + normalized.peak.cpm$k4_temp10_2 +normalized.peak.cpm$k4_temp10_3)/3,
        (normalized.peak.cpm$k4_temp14_1 + normalized.peak.cpm$k4_temp14_2 +normalized.peak.cpm$k4_temp14_3)/3,
        (normalized.peak.cpm$k4_temp20_1 + normalized.peak.cpm$k4_temp20_2 +normalized.peak.cpm$k4_temp20_3)/3,
        (normalized.peak.cpm$k4_temp26_1 + normalized.peak.cpm$k4_temp26_2 +normalized.peak.cpm$k4_temp26_3)/3,
        outline=F,col=c("#00FFFF","#6699FF","#CC33FF","#FF0000","blue"),ylab="CPM",names=c("10ºC","14ºC","20ºC","26ºC"),las=2,cex.lab=2)





normalized.peak.names <- paste(paste(normalized.peak.cpm$chr, normalized.peak.cpm$start, sep = "_"), normalized.peak.cpm$end, sep = "_")
normalized.peak.cpm <- as.matrix(normalized.peak.cpm[,4:15])
rownames(normalized.peak.cpm) <- normalized.peak.names
head(normalized.peak.cpm)
dim(normalized.peak.cpm)

library(limma)

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4,4)))
colnames(experimental.design) <- c("k4_temp10","k4_temp14","k4_temp20","k4_temp26")

experimental.design


linear.fit <- lmFit(normalized.peak.cpm, experimental.design)

contrast.matrix <- makeContrasts(k4_temp14-k4_temp10,k4_temp20-k4_temp10,
                                 k4_temp26-k4_temp10,
                                 levels=c("k4_temp10","k4_temp14","k4_temp20","k4_temp26"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)





temp14c.vs.temp10c <- topTable(contrast.results, number=nrow(normalized.peak.cpm),coef=1,sort.by="logFC")
head(temp14c.vs.temp10c)
peak.ids.14c.10c <- rownames(temp14c.vs.temp10c) 
fold.change.14c.10c <- temp14c.vs.temp10c$logFC
q.value.14c.10c <- temp14c.vs.temp10c$adj.P.Val

activated.peaks.temp14c.vs.temp10c <- peak.ids.14c.10c[fold.change.14c.10c > 1 & q.value.14c.10c < 0.05]
repressed.peaks.temp14c.vs.temp10c <- peak.ids.14c.10c[fold.change.14c.10c < - 1 & q.value.14c.10c < 0.05]

length(activated.peaks.temp14c.vs.temp10c)
length(repressed.peaks.temp14c.vs.temp10c)


temp20c.vs.temp10c <- topTable(contrast.results, number=nrow(normalized.peak.cpm),coef=2,sort.by="logFC")
head(temp20c.vs.temp10c)
peak.ids.20c.10c <- rownames(temp20c.vs.temp10c) 
fold.change.20c.10c <- temp20c.vs.temp10c$logFC
q.value.20c.10c <- temp20c.vs.temp10c$adj.P.Val

activated.peaks.temp20c.vs.temp10c <- peak.ids.20c.10c[fold.change.20c.10c > 1 & q.value.20c.10c < 0.05]
repressed.peaks.temp20c.vs.temp10c <- peak.ids.20c.10c[fold.change.20c.10c < - 1 & q.value.20c.10c < 0.05]

length(activated.peaks.temp20c.vs.temp10c)
length(repressed.peaks.temp20c.vs.temp10c)



temp26c.vs.temp10c <- topTable(contrast.results, number=nrow(normalized.peak.cpm),coef=3,sort.by="logFC")
head(temp26c.vs.temp10c)
peak.ids.26c.10c <- rownames(temp26c.vs.temp10c) 
fold.change.26c.10c <- temp26c.vs.temp10c$logFC
q.value.26c.10c <- temp26c.vs.temp10c$adj.P.Val


activated.peaks.temp26c.vs.temp10c <- peak.ids.26c.10c[fold.change.26c.10c > 1 & q.value.26c.10c < 0.05]
repressed.peaks.temp26c.vs.temp10c <- peak.ids.26c.10c[fold.change.26c.10c < - 1 & q.value.26c.10c < 0.05]

length(activated.peaks.temp26c.vs.temp10c)
length(repressed.peaks.temp26c.vs.temp10c)

boxplot(rowMeans(2^normalized.peak.cpm[activated.peaks.temp26c.vs.temp10c,paste("k4_temp10",1:3,sep="_")]),
        rowMeans(2^normalized.peak.cpm[activated.peaks.temp26c.vs.temp10c,paste("k4_temp14",1:3,sep="_")]),
        rowMeans(2^normalized.peak.cpm[activated.peaks.temp26c.vs.temp10c,paste("k4_temp20",1:3,sep="_")]),
        rowMeans(2^normalized.peak.cpm[activated.peaks.temp26c.vs.temp10c,paste("k4_temp26",1:3,sep="_")]),outline=F,col=c("#00FFFF","#6699FF", "#CC33FF","#FF0000"),ylab="CPM",names=c("10ºC","14ºC","20ºC","26ºC"),las=2,cex.lab=1)


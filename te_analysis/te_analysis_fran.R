ltr <- read.table(file="ltr_te_elements_out.txt",header=F,as.is=T)
head(ltr)
nrow(ltr)

dna.transposons <- read.table(file = "dna_transposons_te_elements_out.txt",
                              header=F, 
                              as.is=T)
head(dna.transposons)
nrow(dna.transposons)

write.table(x = dna.transposons,file = "dna_transposons_unsorted.bed",sep = "\t",row.names = F,col.names = F)


ltr.transposons <- read.table(file = "ltr_te_elements_out.txt",
                              header=F, 
                              as.is=T)
head(ltr.transposons)
nrow(ltr.transposons)

write.table(x = ltr.transposons,file = "ltr_transposons_unsorted.bed",sep = "\t",row.names = F,col.names = F)



unknown.transposons <- read.table(file = "unknown_te_elements_out.txt",
                              header=F, 
                              as.is=T)
head(unknown.transposons)
nrow(unknown.transposons)

write.table(x = unknown.transposons,file = "unknown_transposons_unsorted.bed",sep = "\t",row.names = F,col.names = F)

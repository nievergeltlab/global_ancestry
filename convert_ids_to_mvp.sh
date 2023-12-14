#Datasets like the MVP do not have standard marker names, here is how ot update the SNpweights file:

head -n 5  gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel >  gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel.head
tail -n+6 gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel >  gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel.tail

R

library(data.table)
conversion <- fread('MVP_postqc.snp_name_translation.txt',data.table=F)
names(conversion) <- c("AFFY","SNP")

#Do snplist
snplist <- fread('gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snplist',data.table=F,header=F)
names(snplist) <- c("SNP")
snplist_alt <- subset(merge(snplist,conversion,by="SNP"),select="AFFY")

#Do SNPweights
snpweights <- fread('gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel.tail',data.table=F,header=F)
names(snpweights)[1] <- c("SNP")
snpweights$order <- 1:nrow(snpweights)
snpweights2 <- merge(conversion,snpweights,by="SNP")

snpweights2$SNP <-snpweights2$AFFY
snpweights2$AFFY <- NULL
snpweights2 <- snpweights2[order(snpweights2$order),]
snpweights2$order <- NULL

write.table(snpweights2,file="gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel.tail.mvp",quote=F,row.names=F,col.names=F)
write.table(snplist_alt,file="gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snplist.mvp",row.names=F,quote=F,col.names=F)
q() 
n

#Make MVP ready file
cat  gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel.head gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel.tail.mvp  > gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel.mvp




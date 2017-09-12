#Rename duplicate rs-ids in a .bim file
#Sometimes you get .bim files where rs-ids are not unique. This renames the duplicates. Preference to not rename is given in order of first marker listed.
#V1 - Oct 28, 2016 
args <- commandArgs(trailingOnly = TRUE)
qced_bim <- args[1]


print("Reading .bim file")
 dat1 <- read.table(qced_bim,stringsAsFactors=F,header=F) #Load the QCed SNP list
 names(dat1) <- c("CHR","SNP","cM","BP","A1","A2")


 dat1$SNP_rename <- dat1$SNP #Make a new entry of what the SNP should be called

print("Renaming duplicates")
 #rename the duplicated markers, run three times to make sure multiple duplications are removed
 try(dat1[duplicated(dat1$SNP_rename),]$SNP_rename <- paste(dat1[duplicated(dat1$SNP_rename),]$SNP_rename,"_duplicate",sep=""))
 try(dat1[duplicated(dat1$SNP_rename),]$SNP_rename <- paste(dat1[duplicated(dat1$SNP_rename),]$SNP_rename,"_duplicate",sep=""))
 try(dat1[duplicated(dat1$SNP_rename),]$SNP_rename <- paste(dat1[duplicated(dat1$SNP_rename),]$SNP_rename,"_duplicate",sep=""))
print("Saving data")
 write.table(subset(dat1,select=c(CHR,SNP_rename,cM,BP,A1,A2)),paste(qced_bim,"_rename",sep=""), row.names=F,quote=F,col.names=F)

args <- commandArgs(trailingOnly = TRUE)
qced_bim <- args[1]
missingness <- args[2]
rsid_update <- args[3]

print("Reading .bim file")
 dat1a <- read.table(qced_bim,stringsAsFactors=F,header=F) #Load the QCed SNP list
 names(dat1a)[2] <- "SNP"
print("Reading missingness file")
 dat1b <- read.table(missingness,stringsAsFactors=F,header=T) #Load the missingness per SNP

print("Reading rs-id file")
 refs <- read.table(rsid_update, stringsAsFactors=F,header=T) #Load the b138 SNP names (ones renamed to '.' have been excluded from this list. SNPs with multiple entries take the first entry only)
 names(refs)[1] <- "SNP"
 names(refs)[2] <- "SNP_fix"

print("Merging datasets")
 dat1c <- merge(dat1a,dat1b,by="SNP",all.x=TRUE) #Pool SNP and missingness data. all.x=true is a precauation, the files should be the same dimension and ordering
 dat1 <- merge(dat1c,refs,by="SNP",all.x=TRUE) #Load the updated SNP names
 dat1[is.na(dat1$SNP_fix), ]$SNP_fix <- dat1[is.na(dat1$SNP_fix), ]$SNP #For custom SNPs, make the SNP_fix column just the SNP column
print("Identifying duplicates")
 dat1 <- dat1[order(dat1$N_MISS),] #Order by missingness from least to greatest (later order entries are noted as duplicates)
 dat1$SNP_rename <- dat1$SNP_fix #Make a new entry of what the SNP should be called

print("Renaming duplicates")
 #rename the duplicated markers, run three times to make sure multiple duplications are removed
 try(dat1[duplicated(dat1$SNP_rename),]$SNP_rename <- paste(dat1[duplicated(dat1$SNP_rename),]$SNP_rename,"_duplicate",sep=""))
 try(dat1[duplicated(dat1$SNP_rename),]$SNP_rename <- paste(dat1[duplicated(dat1$SNP_rename),]$SNP_rename,"_duplicate",sep=""))
 try(dat1[duplicated(dat1$SNP_rename),]$SNP_rename <- paste(dat1[duplicated(dat1$SNP_rename),]$SNP_rename,"_duplicate",sep=""))
print("Saving data")
 write.table(subset(dat1,select=c(SNP,SNP_rename)),paste(qced_bim,"_allele_update",sep=""), row.names=F,quote=F)

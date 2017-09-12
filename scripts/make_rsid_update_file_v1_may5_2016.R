#Update Illumina rs-id conversion sheet to be useful for PLINK analysis

#V1 - May 5, 2016 - Initial release

args <- commandArgs(trailingOnly = TRUE)
rsid_illumina <- args[1]

dat1 <- read.table(rsid_illumina,stringsAsFactors=F,header=T)
names(dat1)[2] <- "SNP"

dat1 <- subset(dat1,SNP != ".") #From this list, remove the markers that rename to "." . The markers may not exist anymore but it's not helpful at all to rename them to "."




    unlist_split <- function(x, ...)
	{
		toret <- unlist(strsplit(x, ...))[1]
		return(t(toret))
	}
dat1$SNP2 <- NA

#If a SNP has multiple entries, which will be comma delimited, show the first entry only
ns <- t(t(sapply(dat1$SNP,unlist_split,split=",")))

dat1$SNP2 <- t(t(sapply(dat1$SNP,unlist_split,split=",")))

nm <- which(dat1$SNP2 != dat1$SNP)	

print(paste(length(nm), 'SNPs have multiple names. Picking first name...'))

#This is a list of the markers that can be renamed, with their new names
write.table(dat1[,c("Name","SNP2")],paste(rsid_illumina,"nodot_first",sep="_"), row.names=F,quote=F)

#rename the duplicated markers. Run this 3 times over to get rid of double, and even triple duplicated markers
try(dat1[duplicated(dat1$SNP2),2] <- paste(dat1[duplicated(dat1$SNP2),2],"_duplicate",sep=""))
try(dat1[duplicated(dat1$SNP2),2] <- paste(dat1[duplicated(dat1$SNP2),2],"_duplicate",sep=""))
try(dat1[duplicated(dat1$SNP2),2] <- paste(dat1[duplicated(dat1$SNP2),2],"_duplicate",sep=""))

#This is a list of all markers that can be renamed, where duplicate markers have been given unique names (this prevents PLINK errors); Note, automatically made in make_allele_update_file.R script. Exists just for user reference!
write.table(dat1[,c("Name","SNP2")],paste(rsid_illumina,"use",sep="_"), row.names=F,quote=F)



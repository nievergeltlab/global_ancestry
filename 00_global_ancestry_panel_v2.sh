### Global ancestry determination

## HGDP +KGP data curation 
# Download all HGDP+KGP data from GNomad ( https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg ) 
for chr in {1..22}
do
 wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr"$chr".vcf.bgz
 wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr"$chr".vcf.bgz.tbi
done

# newrefs.subjects is the list of subjects in this data who belong to a reference population 
# nefrefs.subjects.plink is the same list, in PLINK format (FID and IID)
# newref.ind_fixed is the same list,  but formatted for SNPweights - has ID, SEX, and putative ancestral population
# newref.header is a header with FID IID SEX POPULATION

# This is the original curated list from Caroline for the AIMs paper, minus whatever subjects were not sequenced (paper was based on the original HGDP array data)

# User: You must make a SNP list of all markers in your test dataset. 
 # Example: Make a SNP list from a .bim file
 #For MVP, see script of getting post qc markers. MVP_postqc.snplist

# Subset the HGDP+KGP VCF file to only reference subjects, and only the test data SNP list 
# Then convert data into PLINK format
for chr in {1..22}  # {1..21} #22 14 remains!
do
 #bcftools view -S newref.subjects  --include ID==@MVP_postqc.snplist -m2 -M2  -Oz gnomad.genomes.v3.1.2.hgdp_tgp.chr"$chr".vcf.bgz > gsa/gnomad.genomes.v3.1.2.hgdp_tgp.chr"$chr".FILTERED.USE.vcf.bgz
 plink --vcf gsa/gnomad.genomes.v3.1.2.hgdp_tgp.chr"$chr".FILTERED.USE.vcf.bgz --keep newref.subjects.plink --make-bed --out gsa/gnomad.genomes.v3.1.2.hgdp_tgp.chr"$chr".FILTERED.USE.vcf.bgz
done


# Attempt to merge all datasets. 
# It will very likely give an error about duplicate markers
 ls gsa/* | grep bed | sed 's/.bed//g' | awk '{print $1".bed",$1".bim",$1".fam"}' | grep -v allchr | grep -v USE2 > gsa/gsa.mergelist

 plink --merge-list gsa/gsa.mergelist --make-bed --out gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE

# Remove duplicate markers and attempt to merge again

for chr in {1..22}
do
 plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.chr"$chr".FILTERED.USE.vcf.bgz --set-missing-var-ids unk@_# --exclude gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE-merge.missnp --make-bed --out gsa/gnomad.genomes.v3.1.2.hgdp_tgp.chr"$chr".FILTERED.USE2.vcf.bgz
done

ls gsa/* | grep bed | sed 's/.bed//g' | awk '{print $1".bed",$1".bim",$1".fam"}' | grep -v allchr | grep  USE2 > gsa/gsaF.mergelist

# Notice that in the merge step, the missing phenotype code is 1. This is because eigenstrat will remove data with phenotype values = 0
 plink --merge-list gsa/gsaF.mergelist --make-bed --keep newref.subjects.plink --output-missing-phenotype 1 --out gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA


#Get ambiguous alleles
grep -P "A\tT" gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA.bim > ambiguous_snps.txt
grep -P "T\tA" gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA.bim >> ambiguous_snps.txt
grep -P "C\tG" gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA.bim >> ambiguous_snps.txt
grep -P "G\tC" gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA.bim >> ambiguous_snps.txt

#Need markers with maf > 1% in continental populations
sed 's/:/ /g' newref.ind_fixed | awk '{if ($4 == "EUR")print $1,$2}'  > newref.ind_fixed.eur
sed 's/:/ /g' newref.ind_fixed | awk '{if ($4 == "AFR")print $1,$2}'  > newref.ind_fixed.afr
sed 's/:/ /g' newref.ind_fixed | awk '{if ($4 == "SAS")print $1,$2}'  > newref.ind_fixed.sas
sed 's/:/ /g' newref.ind_fixed | awk '{if ($4 == "EAS")print $1,$2}'  > newref.ind_fixed.eas
sed 's/:/ /g' newref.ind_fixed | awk '{if ($4 == "AMR")print $1,$2}'  > newref.ind_fixed.amr

#List only rare markers for each group (MAF < 1%
plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA --keep newref.ind_fixed.amr --max-maf 0.01 --write-snplist --out  gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER.amr
plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA --keep newref.ind_fixed.afr --max-maf 0.01 --write-snplist --out  gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER.afr
plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA --keep newref.ind_fixed.sas --max-maf 0.01 --write-snplist --out  gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER.sas
plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA --keep newref.ind_fixed.eas --max-maf 0.01 --write-snplist --out  gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER.eas
plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA --keep newref.ind_fixed.eur --max-maf 0.01 --write-snplist --out  gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER.eur

cat gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER.*.snplist > remove_maf01.snplist

cat ambiguous_snps.txt remove_maf01.snplist > exclusions.snplist 

#filter to only decently genotyped markers, that are not on the excluded marker list
plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USEA --geno 0.05 --exclude exclusions.snplist --make-bed --out  gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER

#LD prune markers
plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER --indep-pairwise 1000 50 0.2
plink --bfile gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE_PREFILTER --extract plink.prune.in --make-bed --out gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE

## SNPweights panel curation

# Set path to convertf/smartpca programs
 convertf_path=/home/genetics/EIG-master/bin/convertf
 smartpca_path=/home/genetics/EIG-master/bin/smartpca

#Note: To run SNPweights, you need a version of smartpca that gives the trace in the log file (5.0.2 and earlier OR recompile)
 
#2. Go to directory /src/eigensrc/
#3. Open smartpca.c
#4. Find the string: printf("trace: %9.3f\n", y);
#5. Remove the comment characters (") to make the code an active line.
#In version 7.2.1, the code is located at line #1079. In version 6.1.4, the code is at line #1138.
#The SMARTPCA program still calculates trace, changing the code outputs the value of trace for
#use by SNPweights software.

 
# Make a parameter file for convertf
echo " 
genotypename:    gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.bed
snpname:         gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.bim # or example.map, either works 
indivname:       gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.fam # or example.ped, either works
outputformat:    EIGENSTRAT
genotypeoutname: gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.geno
snpoutname:      gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snp
indivoutname:    gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.ind
familynames:     YES
outputgroup:	YES
" > gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.convertf_par

#Run convertf
$convertf_path -p gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.convertf_par

# Make a smartpca parameter file
echo "
genotypename: gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.geno
snpname: gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snp
indivname: gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.ind
evecoutname: gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.evec
evaloutname: gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.eval
altnormstyle: NO
numoutevec: 50
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 6
qtmode: 0
" > gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.eigen_parfile

#Run smartpca
$smartpca_path -p gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.eigen_parfile  > gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.logf

#Make the snpweights .parameter file
echo "
geno: ../gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.geno
snp:  ../gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snp
ind:  ../newref.ind_fixed
evec: ../gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.evec
eval: ../gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.eval
log:  ../gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.logf
snpwtoutput: ../gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel
" > gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightpar

#Calculate SNPweights for the reference data

#Highly recommend creating a conda environment for SNPweights because of the packages required
 conda activate snpweights
 cd SNPweights2.1
 python2.7 calc_snpwt.py --par ../gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightpar
 cd ..

#Given eigenstrat reference data outputs, reformat this data for plotting purposes
 #Print only relevant columns
 #Fix the header
 #Remove the word Case from output
 #Sort data 
 #Add header and populations origins to data
 
 #adjust EVEC file for plotting
 awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.evec | sed 's/#eigvals:/FID IID /g' |  sed 's/:/ /g' | sed 's/Case//g' | sort -g -k 1,1b | join  <(sed 's/:/ /g' newref.ind_fixed | cat newref.header - | LC_ALL=C sort -k1b,1)  - > gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.evec_fixed
 

R
 #WR colorlist
	colorlist <- rep(rgb(0,0,0,1), 8)
	colorlist[1] <- rgb(242,24,3,255, maxColorValue = 255) #1 are africans
	colorlist[2] <- rgb(0,6,253,255, maxColorValue = 255) #2 europeans
	colorlist[3] <- rgb(142,0,117,255, maxColorValue = 255) #3 han chinese
	colorlist[4] <- rgb(178,178,55,255, maxColorValue = 255) #4 sw asians
	colorlist[5] <- rgb(68,26,80,255, maxColorValue = 255)   #5 native americans
	colorlist[6] <- rgb(0,153,229,255, maxColorValue = 255) #6 oceanians
	colorlist[7] <- rgb(255,0,249,255, maxColorValue = 255) #7 cs asian
	colorlist[8] <- rgb(122,124,166,255, maxColorValue = 255) #8 Admixed

  	  
	refdat_pos <- read.table('gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.evec_fixed', header=T)
	names(refdat_pos)[6:11] <- c("PC1","PC2","PC3","PC4","PC5","PC6")

	refdat_pos$color <- NA
	refdat_pos[which(refdat_pos$WR == "AFR"),]$color <- colorlist[1]
	refdat_pos[which(refdat_pos$WR == "EUR"),]$color <- colorlist[2]
	refdat_pos[which(refdat_pos$WR == "EAS"),]$color <- colorlist[3]
	refdat_pos[which(refdat_pos$WR == "SAS"),]$color <- colorlist[7]
	refdat_pos[which(refdat_pos$WR == "SWA"),]$color <- colorlist[4]	
	refdat_pos[which(refdat_pos$WR== "AMR"),]$color <- colorlist[5]
	refdat_pos[which(refdat_pos$WR == "OCE"),]$color <- colorlist[6]

 pdf('gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.evec_fixed.pdf',7,7)
 plot(refdat_pos$PC1,refdat_pos$PC2,col=refdat_pos$color,pch=18)

  plot(refdat_pos$PC1,refdat_pos$PC3,col=refdat_pos$color,pch=18)
  plot(refdat_pos$PC1,refdat_pos$PC4,col=refdat_pos$color,pch=18)
  plot(refdat_pos$PC1,refdat_pos$PC5,col=refdat_pos$color,pch=18)
  plot(refdat_pos$PC1,refdat_pos$PC6,col=refdat_pos$color,pch=18)
  
  
 dev.off()
 
 
 library(plyr)
 
 mymedians <- ddply(refdat_pos, ~ WR, colwise(median,.(PC1,PC2,PC3,PC4,PC5)))
 write.csv(mymedians, 'gsa/gnomad_mvp.snpweightrefpanel_clustercenters.csv',row.names=F)

q()
n

### Important output files. Use these files as parameter inputs into the global ancestry panel (ancestry_pipeline.sh) as a custom panel!h


#SNPweights file 
 gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel

#Ancestry specific cluster medians for plotting
 gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snpweightrefpanel_clustercenters.csv

#SNP list of all markers in referenece pnael
#Get SNP list
 awk '{print $2}' gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.bim > gsa/gnomad.genomes.v3.1.2.hgdp_tgp.allchr.FILTERED.USE.snplist


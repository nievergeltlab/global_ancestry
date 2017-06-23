### Find ancestry informative markers, using 1000 Genomes reference data as an example

#Use global ancestry as a phenotype for association analysis in unrelated KGP subjects
#Perform LD based clumping so that in each LD block, the strongest AIMs are selected
#Making sure that no indels or ambiguous markers are used
#Liftover markers (if you care)

#Get subject populations

#Extract only AFR and EA samples from the 1000G VCF data, on non monomorphic markers
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt

awk '{if ($2=="europe") print $1}' 1000g_pops_wrs.txt > 1000g_europeans.subjects
awk '{if ($2=="africa") print $1}' 1000g_pops_wrs.txt > 1000g_africans.subjects

#Make EA/A phenotypes
cat 1000g_europeans.subjects 1000g_africans.subjects > 1000g_europeans_africans.subjects

awk '{if ($2=="europe") print $1,".",1}' 1000g_pops_wrs.txt > 1000g_europeans.pheno
awk '{if ($2=="africa") print $1,".",2}' 1000g_pops_wrs.txt > 1000g_africans.pheno
cat 1000g_europeans.pheno 1000g_africans.pheno > 1000g_europeans_africans.pheno

#Extract each chromosome
for i in  {1..22}
do

 #plinkseq-0.10/pseq /home/genetics/Desktop/psychchip_aims_panel/kgp_extract_snps/vcfs/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz write-ped    --mask indiv=@1000g_europeans_africans.subjects meta.req=EUR_AF:gt:0.02,EUR_AF:lt:.98,AFR_AF:gt:.02,AFR_AF:lt:.98 --out kgp_eurafr_chr$i 
 #./plink --tfile kgp_eurafr_chr$i --pheno 1000g_europeans_africans.pheno --make-bed --out kgp_eurafr_chr$i
 #rm -f "kgp_eurafr_chr$i".tped
 awk '{print $2}' kgp_eurafr_chr$i.bim | uniq -d > kgp_eurafr_chr"$i".triallelic
 awk '{if ((length($5) > 1) || (length($6) > 1)) print $2}' kgp_eurafr_chr$i.bim  > kgp_eurafr_chr"$i".indel
   grep -P "A\tT" kgp_eurafr_chr"$i".bim | awk '{print $2}' > ambiguous_snps_chr"$i".txt
   grep -P "T\tA" kgp_eurafr_chr"$i".bim | awk '{print $2}'>> ambiguous_snps_chr"$i".txt
   grep -P "C\tG" kgp_eurafr_chr"$i".bim | awk '{print $2}'>> ambiguous_snps_chr"$i".txt
   grep -P "G\tC" kgp_eurafr_chr"$i".bim | awk '{print $2}'>> ambiguous_snps_chr"$i".txt
   grep "var" kgp_eurafr_chr"$i".bim | awk '{print $2}'>> ambiguous_snps_chr"$i".txt
 cat kgp_eurafr_chr"$i".triallelic kgp_eurafr_chr"$i".indel ambiguous_snps_chr"$i".txt > remove_snps_chr"$i".txt

 ./plink --bfile kgp_eurafr_chr$i --exclude remove_snps_chr"$i".txt  --make-bed --out fkgp_eurafr_chr$i
done


#Attempt to merge the data
ls fkgp_eurafr_chr*.bed | sed 's/.bed//g' | awk '{print $1".bed",$1".bim",$1".fam"}' > kgp_eurafr.mergelist
./plink --merge-list kgp_eurafr.mergelist --allow-no-sex --make-bed --out kgp_eurafr_allchr #The merge step is what identifies tri allelics, because tri allelic markers get listed twice


#Run relatedness on EA and AA samples to filter out relateds
 ./plink --bfile kgp_eurafr_allchr --maf 0.05 --indep-pairwise 50 5 0.2 --out "kgp_eurafr_allchr"_ibd
 ./plink --bfile kgp_eurafr_allchr --mind 0.02 --extract  "kgp_eurafr_allchr"_ibd.prune.in --genome --min 0.2 --out  kgp_eurafr_allchr_ibd
 awk '{if(NR ==1) print "FID","IID"; else print $3,$4}' kgp_eurafr_allchr_ibd.genome > "kgp_eurafr_allchr"_ibd.remove

#Run assoc analysis on race
 ./plink --bfile  kgp_eurafr_allchr --remove "kgp_eurafr_allchr"_ibd.remove --logistic  --allow-no-sex --out kgp_eurafr_allchr_norel_assn


 ./plink --bfile kgp_eurafr_allchr --clump kgp_eurafr_allchr_norel_assn.assoc.logistic --allow-no-sex --remove "kgp_eurafr_allchr"_ibd.remove --clump-p1 5e-7 --clump-p2 5e-7 --clump-r2 .1 --clump-kb 100 --out kgp_eurafr_clump

 awk 'NR>1{print $3}' kgp_eurafr_clump.clumped | grep -v '^ *$'  > kgp_eurafr_clump.snplist

#Take random 100k out of this
 sort --random-sort -k1 kgp_eurafr_clump.snplist | head -n 100000 > kgp_eurafr_clump.snplist.random
 
 ./plink --bfile kgp_eurafr_allchr --remove "kgp_eurafr_allchr"_ibd.remove --extract kgp_eurafr_clump.snplist.random --freq --allow-no-sex  --assoc --out 100k_afs
 ./plink --bfile kgp_eurafr_allchr --remove "kgp_eurafr_allchr"_ibd.remove --extract kgp_eurafr_clump.snplist.random --freq --allow-no-sex  --out 100k_afs




R
d1 <- read.table('kgp_eurafr_clump.snplist',header=F,stringsAsFactors=F,nr=250000)
indel <- which(nchar(d1[,5]) > 1 | nchar(d1[,6]) > 1)
d2<- d1[-indel,]
set.seed(17)
keep <- sample(1:dim(d2)[1],100000,replace=F)
write.table(d2[keep,2],"100k_afdifmarkers.snplist",quote=F,row.names=F)
q("no")



#Concatenate all clumps, filter genotype data to them
cat kgp_eurafr_clump_chr*.clumped | awk '{print $3}' | grep -v SNP > kgp_eurafr_clump_allchr.clumped


 ./plink --bfile kgp_eurafr_allchr  --extract kgp_eurafr_clump_allchr.clumped --make-bed --out kgp_eurafr_allchr_200kclump

#Check for liftover - #Every single one is unchaged evidently 
python liftover.py  kgp_eurafr_allchr_200kclump.map > kgp_eurafr_allchr_200kclump.lift

#Read in the list and remove the indels, filter to 100k markers
R
d1 <- read.table('kgp_eurafr_allchr_200kclump.bim',header=F,stringsAsFactors=F,nr=250000)
indel <- which(nchar(d1[,5]) > 1 | nchar(d1[,6]) > 1)
d2<- d1[-indel,]
set.seed(17)
keep <- sample(1:dim(d2)[1],100000,replace=F)
write.table(d2[keep,2],"100k_afdifmarkers.snplist",quote=F,row.names=F)

#Get basic summary data on these markers...
./plink --bfile kgp_eurafr_allchr  --extract 100k_afdifmarkers.snplist --freq --allow-no-sex  --assoc --out 100k_afs


#With r2=.2, approx 6% of the 6.5 million variants are clumped... want even less coverage

#Update to b144
wget http://genome.sph.umich.edu/wiki/LiftRsNumber.py
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/b144_SNPChrPosOnRef_107.bcp.gz
python liftover.py kgp_eurafr_allchr.map 


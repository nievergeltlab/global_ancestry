###Ricopili Modified for Ancestry (MANC) master script

###Initial configuration steps

#User: Unzip the data, put contents into a working directory

#User:Set path to working directory (i.e. this is where ancestry_pipeline.sh is stored)
 WORKING_DIR="/home/maihofer/MRSA/scripts/bbanctest"
 
##Call into WD
 cd $WORKING_DIR
 

#User: Give the location of the PLINK 2 binary
 plink_location="$WORKING_DIR"/plink

#User: Write the name of the PLINK bed/bim/fam 
 bfile=mrs_gwas1_vPGC_use

#User: Give the folder where this PLINK binary is stored 
 bfile_directory=/home/maihofer/MRSA/starting_data

#User: Location of SNPweights (Download from http://www.hsph.harvard.edu/alkes-price/software/)
 snpweights_path="$WORKING_DIR"/SNPweights2.1/inferancestry.py
 
#User: Location of convertf tool from EIGENSOFT (Download from http://www.hsph.harvard.edu/alkes-price/software/)
 eigensoft_loc="$WORKING_DIR"/EIG5.0.2/bin/convertf
 

#User: Name of panel chosen ('Default': works for most illumina. AffyBB: Designed for Affy biobank chip. GSA: Designed for Illumina GSA)
 panel=Default
 
#This will supply: 
#Name of the list of ancestry panel SNP rsids (packaged with script, be sure that this is in WORKING_DIR)
#Name of ancestry panel (packaged with script, be sure that this is in WORKING_DIR)
#Name of ancestry panel cluster centers (packaged with script, be sure that this is in WORKING_DIR)
 
 if [ $panel == "Default" ]
 then
  snpweights_snplist=hgdp_kgp_merged_v3_jointsample_v4_k6.snplist
  snpweightfile_path=hgdp_kgp_merged_v3_jointsample_v4_k6.snpweightrefpanel
  snpweight_clustercenters=hgdp_kgp_merged_v3_jointsample_v4_k6.snpweightrefpanel_clustercenters.csv
 fi
 
  if [ $panel == "GSA" ]
 then
  snpweights_snplist=hgdp_kgp_merged_v3_jointsample_v4_k6.snplist
  snpweightfile_path=hgdp_kgp_merged_gsa_v3_jointsample_v4_k6.snpweightrefpanel
  snpweight_clustercenters=hgdp_kgp_merged_gsa_v3_jointsample_v4_k6_forsnpweights.snpweightrefpanel_clustercenters.csv
 fi
 
 if [ $panel == "AffyBB" ]
 then
  snpweights_snplist=hgdp_kgp_merged_v3_jointsample_v4_k6_bb.snplist
  snpweightfile_path=hgdp_kgp_merged_v3_jointsample_v4_k6_bb.snpweightrefpanel
  snpweight_clustercenters=hgdp_kgp_merged_v3_jointsample_v4_k6_bb.snpweightrefpanel_clustercenters.csv
 fi
 
#Give path of Manufacturer SNP ID to rs-id conversion table.
#Illumina/Affy provide these 
#UK Biobank Array snplist already provided

#For Illumina panels: Go to the kit support page for the array, click downloads, then click ArrayNameHere support files
#the file should be under the name of "ArrayNameHere Loci Name to rsID Conversion File"

 illumina_snplist=Multi-EthnicGlobal_B1_b144_rsids.txt # ukb_affy_to_snp.txt
#multi-ethnic-global-8-b1-rsids.zip
###Ancestry determination steps

#Use Illumina supplied list of SNPs that can be renamed
 Rscript --vanilla scripts/make_rsid_update_file_v1_may5_2016.R "$illumina_snplist"

 rsidfile="$WORKING_DIR"/"$illumina_snplist"_nodot_first

## Call subject ancestries (1000 subjects are processed at a time enhance built in parallelization)
#Note: This is programmed to run on a local node. Generally this is not a problem even for hundreds of thousands of subjects.
#If you MUST run this as a job script, then edit the call_ancestry_v2_may6_2016.sh file,
#replacing the final 'bash' command with your job submission command.
 chmod u+rwx scripts/call_ancestry_v2_may6_2016.sh
 bash scripts/call_ancestry_v2_may6_2016.sh $plink_location $bfile $rsidfile $bfile_directory $snpweights_snplist $snpweights_path $snpweightfile_path $eigensoft_loc

#Combine the SNPweights ancestry calls
 cat temporary_files/"$bfile"_anc_*.predpc_oneweek  > ancestry/"$bfile".predpc_oneweek

#Call into the ancestry folder, classify subject ancestry, produce an ancestry PC plot
 cd "$WORKING_DIR"/ancestry 
 Rscript --vanilla "$WORKING_DIR"/scripts/ancestry_plots_v4_mar1_2017.R "$bfile".predpc_oneweek ../"$snpweight_clustercenters" $panel

 



#!/bin/bash
#Call ancestry using SNPweights and reference panel

#V2 - May 6, 2016 - Get rid of cm positions (porder check fails if cm order does not match phys dist order)
#V1 - May 5, 2016 - Initial version


#Syntax:
#call_ancestry plink_location bfile rsidfile bfile_directory snpweights_snplist snpweights_path snpweightfile_path eigensoft_loc

echo "Plink Binary Location: $1"
echo "Bfile Name: $2"
echo "Bfile location: $4"
echo "SNP rsid name location: $3"
echo "SNP list location: $5"
echo "SNPweights inferanc.py location: $6"
echo "SNPweights reference panel location: $7"
echo "convertf location: $8"


#Make a temporary_files directory if not existing
if [ ! -d "temporary_files" ]; then
  mkdir temporary_files
fi


#Run PLINK. Update names of markers to rs-ids. Keep only relevant variants. Get rid of monomorphic and non-autosomal variants. 
echo "Updating marker names"
$1 --bfile "$4/"$2 --update-name $3 --alleleACGT --autosome --maf 0.00000001 --make-bed --out temporary_files/"$2"_anc

#Get rid of cm positions
mv temporary_files/"$2"_anc.bim temporary_files/"$2"_anc.bim.bk
awk '{print $1,$2,"0",$4,$5,$6}' temporary_files/"$2"_anc.bim.bk > temporary_files/"$2"_anc.bim

#Split the .fam files into segments for ancestry calling
echo "Splitting .fam file"
awk '{print $1,$2}' temporary_files/"$2"_anc.fam | split  -l 1000 -d - temporary_files/"$2"_anc_

echo "Converting to Eigenstrart format"
#Convert all data to eigenstrat format
for i in {00..99}
do 
     if [ -e temporary_files/"$2"_anc_"$i" ] 
      then	
           #Subset to only the N subject segment.  Extract relevant AIMs here, otherwise duplicate marker error occurs
           $1  --bfile temporary_files/"$2"_anc --keep temporary_files/"$2"_anc_"$i" --output-missing-phenotype 2 --extract $5 --make-bed --out temporary_files/"$2"_anc_"$i"

           #Make a parameter file for convertf to get these subjects into eigenstrat format
            echo " 
            genotypename:    temporary_files/"$2"_anc_"$i".bed
            snpname:         temporary_files/"$2"_anc_"$i".bim # or example.map, either works 
            indivname:       temporary_files/"$2"_anc_"$i".fam # or example.ped, either works
            outputformat:    EIGENSTRAT
            genotypeoutname: temporary_files/"$2"_anc_"$i".geno
            snpoutname:      temporary_files/"$2"_anc_"$i".snp
            indivoutname:    temporary_files/"$2"_anc_"$i".ind
            familynames:     YES
            outputgroup:	YES
            #pordercheck: NO
            " > temporary_files/"$2"_anc_"$i".convertf_par
          
            $8 -p temporary_files/"$2"_anc_"$i".convertf_par
        fi
done



echo "Writing jobs"
#Construct the job code. Assume 16 core nodes by default

cwd=$(echo pwd)

job=1
nodesize=16
for i in $(seq -w 00 $nodesize 99) #{00..100..$nodesize}
do 
    if [ -e temporary_files/"$2"_anc_"$i" ]
    then
     echo "

        #PBS -N ancestry_determination
        #PBS -l nodes=1:ppn=16
        #PBS -l walltime=2:00:00
        #PBS -j oe

       cd $PWD
          " > temporary_files/ancjob"$job"_"$2".pbs

    #Now for each element of this job we're going to make parameter files for snpweights and a command to run snpweights
    ix=$(($i+$nodesize))

        for j in $(seq -w $i 1 $ix)
        do
          if [ -e temporary_files/"$2"_anc_"$j" ]
           then
          #Make parameter file for SNPweights
echo "
geno:  temporary_files/"$2"_anc_"$j".geno
snp:   temporary_files/"$2"_anc_"$j".snp
ind:   temporary_files/"$2"_anc_"$j".ind
snpwt: $7
predpcoutput: temporary_files/"$2"_anc_"$j".predpc_oneweek
" > temporary_files/"$2"_anc_"$j".snpweights_parfile 

            #Run snpweights to obtain ancestry estimates for test subjects
            echo "$6 --par temporary_files/"$2"_anc_"$j".snpweights_parfile &" >> temporary_files/ancjob"$job"_"$2".pbs
           fi
        done
     echo "wait" >> temporary_files/ancjob"$job"_"$2".pbs #Add "wait" to the end so that all jobs finish prior to the script closing
     job=$((job+1))

    fi
done

totjobs=$(($job - 1))
for i in  $(seq 1 1 $totjobs)
 do 
 echo "Submitting job"
 bash temporary_files/ancjob"$i"_"$2".pbs
done


#!/bin/bash
## number of cores
#SBATCH -n 8
#SBATCH -o macs2-test2.out


#----
# conda activate ATAC-seq
#----

#-----
#-- PATHS
#-----

#--- Input files path --#
FASTQ_DIR="$1"
echo $FASTQ_DIR

REF_DIR="$2"
GENOME_DIR=${REF_DIR}/GRCm38_bwa_index
BLACKLIST_DIR=${REF_DIR}

#FASTQ_DIR=/rugpfs/fs0/tavz_lab/scratch/schhabria/20240104_Mira_ATACseq/fastq/test
#GENOME_DIR=/rugpfs/fs0/tavz_lab/scratch/schhabria/ref_files/mouse_ref/GRCm38-mm10/GRCm38_bwa_index

#--- Results paths --#
RESULTS_DIR="$3"
#RESULTS_DIR=/rugpfs/fs0/tavz_lab/scratch/schhabria/20240104_Mira_ATACseq/outputs-all

#--- Prepreprocessing paths
FASTQC_DIR=$RESULTS_DIR/fastQC
NGmerge_TRIM_FASTQ_DIR=$RESULTS_DIR/trim_ngmerge_fastq

FASTQC_TRIM_NGMERGE_DIR=$RESULTS_DIR/trim_ngmerge_fastQC
FASTQC_TRIM_MIRA_DIR=$RESULTS_DIR/trim_fastQC_Mira

ALIGN_DIR=$RESULTS_DIR/align

#--- Filtering paths
FILTER_DIR=$RESULTS_DIR/filter
STATS_DIR=$FILTER_DIR/stats_reads
MITO_DIR=$FILTER_DIR/mito_rm_bams
DEDUP_DIR=$FILTER_DIR/dedup_bams
ATACQC_DIR=$FILTER_DIR/ATAC_QC
FINAL_FILTER_DIR=$FILTER_DIR/final_filter_bams

#--- Downstream paths
BW_DIR=$RESULTS_DIR/coverage_files
PEAKS_DIR=$RESULTS_DIR/peak_files
MACS2_DIR=$PEAKS_DIR/macs2_peaks
GENRICH_DIR=$PEAKS_DIR/genrich_peaks


#-----
# Parameters
#-----

species=$4

#----**
# *---------- PART 1: PRE_PROCESSING --------*
# ** FastQC
# ** Trimming adapters
# ** FastQC on trimmed fastq
# ** Alignment and sorting
#----**

# ----
# Step 1: Run fastQC on the fastq files 
# ----
# Input: fastq.gz
# Output: fastqc.html
#----- 

cd $FASTQ_DIR
for fastq in *_S*_R*_001.fastq.gz; 
   do  
    samplename=${fastq%.fastq.gz}
    echo "Processing sample ${fastq}"
    echo "Sample name ${samplename}"
    mkdir -p $FASTQC_DIR
    fastqc -o $FASTQC_DIR $fastq   
done


# # ----
# # Step 2: Adapter Trimming
# # ----
# # Input:  fastq.gz
# # Output: trimmed_fastq.gz
# # Tool : Cutadapt - when the adapter sequence is known
# #        NGmerge  - when the adapter sequence is unknown
# #        -a : adapter removal mode and -z to zip all the files
# #----- 

cd $FASTQ_DIR
for read1 in *_S*_R1_001.fastq.gz; 
   do  
    read2=${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}
    samplename=${read1%_R1_001.fastq.gz}
    
    echo "Processing sample ${read1}"
    echo "Processing sample ${read2}"
    echo "Processing sample ${samplename}" 
    mkdir -p ${NGmerge_TRIM_FASTQ_DIR}
    
    NGmerge -a -1 ${read1} -2 ${read2} -o ${NGmerge_TRIM_FASTQ_DIR}/${samplename}_trimmed -z

done


# # ----
# # Step 3: Repeat the fastq on trimmed files to check the quality
# # ----
# # Input: trimmed_fastq.gz
# # Output: trimmed_fastq.html
# #----- 

cd $NGmerge_TRIM_FASTQ_DIR

for fastq in *_trimmed_*.fastq.gz; 
   do  
    samplename=${fastq%.fastq.gz}
    echo "Processing sample ${fastq}"
    echo "Sample name ${samplename}"

    mkdir -p $FASTQC_TRIM_NGMERGE_DIR
    

    fastqc -o $FASTQC_TRIM_NGMERGE_DIR $fastq   
    
done


# #----
# # Step 4: Align the trimmed fastq files 
# #----
# # Input  : trimmed.fastq.gz
# # Output : sample.bam
# # Tool   : bwa-mem
# #----


cd $NGmerge_TRIM_FASTQ_DIR
for read1 in *_trimmed_1.fastq.gz; 
   do  
    read2=${read1/_trimmed_1.fastq.gz/_trimmed_2.fastq.gz}
    samplename=${read1%_trimmed_1*.fastq.gz}
    
    echo "Processing sample ${read1}"
    echo "Processing sample ${read2}"
    echo "Processing sample ${samplename}"  
    
    mkdir -p $ALIGN_DIR
    bwa mem -t 8 ${GENOME_DIR}/Mus_musculus.GRCm38.dna.primary_assembly.fa \
                 ${read1} ${read2} | samtools sort -o ${ALIGN_DIR}/${samplename}_sorted.bam -

    samtools index ${ALIGN_DIR}/${samplename}_sorted.bam

done

# #----**
# # *---------- PART 2: FILTERATING --------*
# # ** Mitochondria removal
# # ** Mark Duplicates
# # ** Remove duplicated and low reads
# # ** Remove Encode blacklisted  region
# # ** Shift reads
# #----**

# #----
# # Step 1: Remove mitochrondrial reads
# #-----
# # Input  : <sample>_sorted.bam
# # Output : <sample>_sorted.idxstats
# # Tool   : samtools idxstats
# #----

cd $ALIGN_DIR
for bam in *_sorted.bam; 
   do  
    samplename=${bam%_sorted.bam}
    
    echo "Processing bam ${bam}"
    echo "Processing sample ${samplename}"  
    
    mkdir -p ${STATS_DIR}
    mkdir -p ${MITO_DIR}

    #--Get the total percentage
    samtools idxstats ${bam} > ${STATS_DIR}/${samplename}_sorted.idxstats
    grep "chrM" ${STATS_DIR}/${samplename}_sorted.idxstats 

    #--Generate the flagstat report (gives you total number of DNA fragments)
    samtools flagstat ${samplename}_sorted.bam > ${STATS_DIR}/${samplename}_sorted.flagstat
    
    #-- Remove any mitocondrial DNA
    samtools view -h ${bam} | grep -v chrM | samtools sort -O bam -o ${MITO_DIR}/${samplename}_rmChrM.bam -T ${MITO_DIR}

done

# #----
# # Step 2: Mark Duplicates (This step is just marking the duplicates)
# #----
# # Input  : rmChrM.bam
# # Output : marked.bam
# # Tool   : picard
# #          QUIET=true (don't suppress job-summary info on System.err.)
# #          VALIDATION_STRINGENCY=LENIENT ()
# #-------


cd ${MITO_DIR}
for bam in *_rmChrM.bam; 
   do  
    samplename=${bam%_rmChrM.bam}
    
    echo "Processing bam ${bam}"
    echo "Processing sample ${samplename}"  
    
    mkdir -p ${DEDUP_DIR}

    #--Run the picard Markduplicates
    picard MarkDuplicates QUIET=true \
           INPUT=${samplename}_rmChrM.bam \
           OUTPUT=${DEDUP_DIR}/${samplename}_marked.bam \
           METRICS_FILE=${DEDUP_DIR}/${samplename}_dup.metrics \
           REMOVE_DUPLICATES=false \
           CREATE_INDEX=true \
           VALIDATION_STRINGENCY=LENIENT \
           TMP_DIR=${DEDUP_DIR}

   #View the % of duplicates
   head -n 8 ${DEDUP_DIR}/${samplename}_dup.metrics | cut -f 7,9 | grep -v ^# | tail -n 2 

done

# #----
# # Step 3: Remove the low quality and duplicates 
# #----
# # Input  : marked.bam
# # Output : filtered.bam
# # Tool   : samtools
# # 2-retain properly mapped pairs; 12-unmapped; 1024-duplicate reads; 512-failed vendor checks >> 1548 - combined all these parameters 
# #-------


cd ${DEDUP_DIR}
for bam in *_marked.bam; 
   do  
    samplename=${bam%_marked.bam}
    
    echo "Processing bam ${bam}"
    echo "Processing sample ${samplename}"  
    

    #-- Filter for PCR duplicates, multi-map reads. low quality based on flags
    samtools view -h -b -f 2 -F 1548 -q 30 ${bam} | samtools sort -o ${samplename}_filtered.bam
    samtools index ${samplename}_filtered.bam  

done

#---
# Step 4: Remove Encode blacklisted 
#---
# Input : 
# Output: 
#----

cd ${DEDUP_DIR}
for bam in *_filtered.bam; 
   do  
    samplename=${bam%_filtered.bam}
    
    echo "Processing bam ${bam}"
    echo "Processing sample ${samplename}"  
    
    #Remove reads within the blacklist regions
    bedtools intersect -nonamecheck -v -abam ${bam} -b ${BLACKLIST_DIR}/mm10-blacklist.v2.bed.gz > ${samplename}.tmp.bam

    #Sort and index the bam file
    samtools sort -O bam -o ${samplename}_blacklist-filtered.bam ${samplename}.tmp.bam
    samtools index ${samplename}_blacklist-filtered.bam

    rm ${samplename}.tmp.bam

done



#----**
# *---------- PART 3: DOWNSTREAM ANALYSIS --------*
# ** Peak Calling with MACS2 and genrich
#----**


# #---
# # Step 1: Peak Calling 
# #---
# # Input  : blacklist-filtered.bam
# # Output : narrowPeak.gz
# # Tool   : MACS2 and Genrich
# # --------


cd ${DEDUP_DIR}
for bam in *blacklist-filtered.bam;
    do  
     samplename=${bam%_blacklist-filtered.bam}
     
     
     mkdir -p ${GENRICH_DIR}
     mkdir -p ${MACS2_DIR}
    
     echo "Processing bam ${bam}"
     echo "Processing sample ${samplename}"  
     
     #-- Method 1: macs2 peakcall
     macs2 randsample -i ${bam} -f BAMPE -p 100 -o ${MACS2_DIR}/${samplename}_macs2.bed

     macs2 callpeak -t ${MACS2_DIR}/${samplename}_macs2.bed \
                    -f BEDPE \
                    --nomodel \
                    --keep-dup all \
                    --nolambda \
                    -g mm \
                    -n ${samplename} --outdir ${MACS2_DIR} 2> ${MACS2_DIR}/macs2.log

     #-- Method 2: genrich peakcall

     samtools view -h ${bam} | samtools sort -n - | samtools view -b - > ${samplename}_nsorted.bam
     Genrich -j -t ${samplename}_nsorted.bam  -o ${GENRICH_DIR}/${samplename}_genrich.narrowPeak
      
    done




#---
# Step 3 --- Annotate the peaks
# -----
# Input  : narrow_peaks.narrowPeak (from macs2)
# Output : annotated_peaks.txt
# Tool   : HOMER
#-----


# cd ${MACS2_DIR}
# #-- Add the Homer annotate_peaks to the path
# PATH=$PATH:/rugpfs/fs0/tavz_lab/scratch/schhabria/Module_Project/eCLIP-ENCODE/HOMER/.//bin/

# for bed in *_narrow_peaks.narrowPeak ; 
#    do  
#     samplename=${bed%narrow_peaks.narrowPeak }
   
#     echo "Processing bed ${bed}"
#     echo "Processing sample ${samplename}"  
    
#     annotatePeaks.pl $bed $species > ${samplename}_annotated.txt

     
# done
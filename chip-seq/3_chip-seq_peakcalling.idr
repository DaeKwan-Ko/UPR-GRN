#------------------------------------------------------------------------------------
# Aim: peaking of cleaned 50-nt single-end ChIP-seq data and obtaining reproducible peaks
# Author: Dae Kwan Ko (dkko@msu.edu)
# Last modified: 01-12-2021
# Usage: shown below
#------------------------------------------------------------------------------------

### merge input files (ENCODE recommends to pool control (input) samples if # reads is less than in ChIP samples)

### peak calling using macs2 in individual replicate using macs2
macs2 callpeak -t ${LINE}_chip.final.bam -c ${INPUT} -f BAM -g 1.19e8 --keep-dup=1 --name ${LINE}.v6 \
--bdg --call-summits --pvalue=1e-2 --nomodel --shift 80 --extsize 160 --outdir /mnt/home/dkko/project_bzip28bzip60/b28.b60_ChIPseq/6_macs2_onlyChr.v6/${LINE}

### peak calling using macs2 in pooled samples using macs2
macs2 callpeak -t ${LINE}_chip_final_pooled.bam -c ${INPUT} -f BAM -g 1.19e8 --keep-dup=1 --name ${LINE}_pooled.v6 \
--bdg --call-summits --pvalue=1e-2 --nomodel --shift 80 --extsize 160 --outdir /mnt/home/dkko/project_bzip28bzip60/b28.b60_ChIPseq/6_macs2_onlyChr.v6/${LINE}_pooled

### idr
idr --samples ${INPUT_1} ${INPUT_2} --peak-list ${POOLED_PEAK_FILE} --input-file-type narrowPeak --output-file ${IDR_OUTPUT} \
--rank signal.value --soft-idr-threshold ${IDR_THRESH} --plot --use-best-multisummit-IDR

### get peaks passing IDR threshold of 5%
awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${IDR_OUTPUT} | sort -k 1,1 -k2,2n | uniq > ${FILT_OUTPUT}

### keep Tm-specific peaks (removed if overlapped >30%) using bedtools
bedtools intersect -v -f 0.30 -a ${Tm_PEAKS} -b ${DMSO_PEAKS} | sort -k 1,1 -k2,2n | uniq > ${OUTPUT}

###  keep Tm-specific peaks (among overlapped peaks, keep if -(log10P-value) of Tm is twice higher than of DMSO)
bedtools intersect -wa -wb -f 0.30 -a ${Tm_PEAKS} -b ${DMSO_PEAKS} | awk 'BEGIN{OFS="\t"} $8>($18*3) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort -k 1,1 -k2,2n | uniq > ${OUTPUT}

### merge filtered peaks in single files
cat ${NON_OVERLAP} ${OVERLAP_3FP} | sort -k 1,1 -k2,2n | uniq > ${OUTPUT}

### merge filtered peak list from different time-points to create union list of peaks for each bZIP28 or bZIP60
cat 28-0.final.Tm_specific.narrowPeak 28-12.final.Tm_specific.narrowPeak 28-24.final.Tm_specific.narrowPeak | sort -k 1,1 -k2,2n | uniq > 28_union.final.Tm_specific.narrowPeak
cat 60-0.final.Tm_specific.narrowPeak 60-12.final.Tm_specific.narrowPeak 60-24.final.Tm_specific.narrowPeak | sort -k 1,1 -k2,2n | uniq > 60_union.final.Tm_specific.narrowPeak

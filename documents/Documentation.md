# Documentation v1.0
ATAC-seq and CUT&Tag-seq data analysis pipeline for Bo Zhang`s lab<br/>Last edit: 10/12/2022<br/>For any question please contact: siyuancheng@wustl.edu

**Outline**
1. Parameters and personalized peakcalling options
2. Output example and annotation
3. Data processing details

<br />

## 1. Parameters and personalized peakcalling options
### 1.1 Parameters
`-h`: help information<br/>
`-d`: data type. *ATAC* for ATAC-seq, *CUTnTag* for CUT&Tag-seq<br/>
`-g`: genome reference. For now the supported genoms are: <mm39/mm10/mm9/hg19/hg38/danRer10/danRer11/rn6/dm6><br/>
`-r`: SE for single-end, PE for paired-end<br/>
`-m`: marker. *narrow* for narrow marker(H3K4me3, H3K27ac, etc), *broad* for broad marker(H3K36me3, etc)<br/>
`-o`: experiment reads file 1 or the SE reads file, must be ended by .fastq or .fastq.gz or .sra<br/>
`-O`: experiment reads file 2 if input PE data, must be ended by .fastq or .fastq.gz<br/>
`-i`: input control read file 1 or the SE reads file (only be used when data type is CUT&Tag-seq), must be ended by .fastq or .fastq.gz<br/>
`-I`: input control read file 2 (only be used when data type is CUT&Tag-seq), must be ended by .fastq or .fastq.gz<br/>
`-s`: (optional) peak calling software. (default: *AIAP* for ATAC-seq data, *MACS2* for CUT&Tag-seq narrow data, *ChromHMM* for CUT&Tag-seq broad data) <br/>
&emsp;&emsp; For ATAC-seq data: <AIAP/MACS2/F-seq/HMMRATAC>; for CUT&Tag-seq data: <AIAP/MACS2/F-seq/HMMRATAC/SEACR/ChromHMM><br/>
`-c`: (optional) specify read length minimum cutoff for methylQA filtering, default: 38<br/>
`-t`: (optional) specify number of threads to use, default: 12<br/>
`--personalize`: (optional) user-defined parameters for peak calling step, default: FALSE<br/>
&emsp;&emsp; if user select using personalized peak calling options (just add `--personalize` on command), user need to define **ALL** the optional parameters for peakcalling software, or user need set parameters as **software-defined default** setting)<br/>

### 1.2 Personalized peak calling options
**1.2.1 MACS2 options:**<br/>
`--keep-dup`: --keep-dup <integer> (Tn5DP default: ATAC--1000, CUT&Tag--1; MACS2-default: 1)<br/>
`--nomodel`: --nomodel (Tn5DP default: ATAC--True, CUT&Tag--False; MACS2-default: False)<br/>
`--shift`: --shift <integer> (Tn5DP default: ATAC--0, CUT&Tag--0; MACS2-default: 0)<br/>
`--extsize`: --extsize <integer> (Tn5DP default: ATAC--150, CUT&Tag--200; MACS2-default: 200)<br/>
`--qvalue`: -q | --qvalue <0-1> (Tn5DP default: ATAC--0.01, CUT&Tag--0.01; MACS2-default: 0.05)<br/>
`--callsummits`: --call-summits (Tn5DP default: False; MACS2-default: False)<br/>
`--broad`: --broad (Tn5DP default: Narrow--False, Broad--True; MACS2-default: False)<br/>
(Note: --call-summits and --broad could NOT be used at the same time)
  
**1.2.2 ChromHMM options:** <br/>
`--binsize`: -b binsize <integer> (Tn5DP default: 200; ChromHMM-default: 200)<br/>
`--seed`: -s seed <integer> (Tn5DP default: 999; ChromHMM-default: 999)<br/>
`--foldthresh`: -f foldthresh <0-1> (Tn5DP default: 0; ChromHMM-default: 0)<br/>
`--poissonthreshold`: -p poissonthreshold <0-1> (Tn5DP default: 0.0001; ChromHMM-default: 0.0001)<br/>
`--zerotransitionpower`: -z zerotransitionpower <integer> (Tn5DP default: 8; ChromHMM-default: 8)<br/>
`--init`: -init <information/random/load> (Tn5DP default: information; ChromHMM-default: information)<br/>
  
**1.2.3 SEACR options** <br/>
`--normalize`: normalize step <norm/non> (Tn5DP default: CUT&Tag--norm)<br/>
`--model`: peakcalling step <relaxed/stringent> (Tn5DP default: CUT&Tag--stringent)<br/>
  
**1.2.4 HMMRATAC options**<br/>
`--upper`: -u <integer> (Tn5DP default: ATAC--35, CUT&Tag--35; HMMRATAC-default: 20)<br/>
`--lower`: -l <integer> (Tn5DP default: ATAC--2, CUT&Tag--2; HMMRATAC-default: 10)<br/>
`--blacklist`: -e <blacklist BED file> (Tn5DP default: ATAC--none, CUT&Tag--none; HMMRATAC-default: none)<br/>
`--score`: --score <max/ave/med/fc/zscore/all> (Tn5DPTn5DP default: ATAC--max, CUT&Tag--max; HMMRATAC-default: max)<br/>
`--bedgraph`: --bedgraph <True/False> (Tn5DP default: ATAC--True, CUT&Tag--True; HMMRATAC-default: False)<br/>

<br />

## 2. Output examples

**2.1 ATAC-seq output**

After running the pipeline, there will be a folder called **Processed_${name}**, all intermediate files and final output files are stored there. And it looks like this:<br/>
<img width="700" alt="cuttag_output" src="https://user-images.githubusercontent.com/81212185/195744041-37ebc3b7-63b8-40f7-89bb-55468284fe63.png">

There should be 9 files in total, and the contents are listed below.
| File name | Content |
| :---      | :---    |
| *QC_data_collection_${name}* | stores all intermediate output files from each step |
| *QC_pipe_processing.log* | record the status of each step |
| *step1.1_${name}_cutadapt_PE.trimlog* | cutadapt trimming report |
| *step2.1_trimed_${name}.bam* | aligned bam file from **BWA MEM**, keep this one only for backup purpose |
| *step3.2_normalized_per_10M_${name}.bigwig* | normalized signal per 10M input single ends for visualization purpose |
| *step3.2_rmbl_${name}.bigwig* | full signal for visualization purpose |
| *step3.3_rmbl_${name}.open.bed* | bed file after quality filtering on aligned bam file and signal shifting, this is the input for **MACS2** |
| *step4.1_peakcalling_${software}_${name}_peaks.narrowPeak* | macs2 output, bed file that records peaks called by the software |

<br />

**2.2 CUT&Tag-seq output**

After running the pipeline, there will be a folder called **Processed_${name}**, all intermediate files and final output files are stored there. And it looks like this:<br/><img width="700" alt="cuttag_output" src="https://user-images.githubusercontent.com/81212185/195747872-7e952dce-81c6-46d4-a35d-9d7df55c5613.png">

There should be 11 files in total, and the contents are listed below.
| File name | Content |
| :---      | :---    |
| *QC_data_collection_${name}* | stores all intermediate output files from each step |
| *QC_IgG_data_collection_${name}* | stores all intermediate output files of IgG input from each step |
| *QC_pipe_processing.log* | record the status of each step |
| *step1.1_${name}_cutadapt_PE.trimlog* | cutadapt trimming report |
| *step2.1_trimed_${name}.bam* | aligned bam file from **BWA MEM**, keep this one only for backup purpose |
| *step3.1_methylQA_${name}.bigwig* | BigWig file generated from **methylQA density** |
| *step3.1_methylQA_${name}.extended.bed* | bed file after quality filtering on aligned bam file and signal shifting, generated by methylQA density, this is the input for **MACS2** |
| *step3.2_rmbl_${name}.bedGraph* | bedGraph file aftering quality filtering and blacklist removal, this is the input for **SEACR** |
| *step3.2_rmbl_${name}.bigwig* | full signal for visualization purpose |
| *step4.1_peakcalling_${software}_${name}.bed* | bed file that records peaks called by the software |

<br />

## 3. Data processing details

### Step1, Pre-alignment

**1.1 Trimming by cutadapt**

Tool: cutadapt v3.5<br/>
Input: fastq file<br/>
Output: trimmed fastq file<br/>
```
 $cutadapt -a $adapter_1 -A $adapter_2 --quality-cutoff=15,10 --minimum-length=25 -o 'step1.1_trimed_'$name'_1.fastq' -p 'step1.1_trimed_'$name'_2.fastq' $read1 $read2 >'step1.1_'$name'_cutadapt_PE.trimlog'
```

**1.2 Fastq file quality assessment**

Tool: FastQC v0.11.7<br/>
Input: trimmed fastq file<br/>
Output: fastqc results
```
[ -f $(ls 'step1.1_trimed_'$name*'.fastq' | head -1) ] && $fastqc -t $threads 'step1.1_trimed_'$name*'.fastq' -o .
```

### Step2, reads alignment and reads distribution

**2.1 BWA MEM aligner**

Tool: bwa v0.7.16a, samtools v1.9<br/>
Input: trimmed fastq file<br/>
Output aligned BAM file
```
$bwa mem -t $threads $bwa_ref 'step1.1_trimed_'$name*'.fastq' | $samtools view -bS - | $samtools sort - -O 'bam' -o 'step2.1_trimed_'$name'.bam' -T temp_aln_input
```

**2.2 reads distributioin**

Tools: samtools v1.9, methylQA v0.2.1<br/>
Input: aligned BAM file
Output: reads distribution count
```
s2.2_distri() {
    echo 'step2.2 distribution...'

    $samtools view -h 'step2.1_trimed_'$name'.bam' >input.sam &&
        awk '$5>0' input.sam | sed '/^@/d' - | cat <(grep '^@' input.sam) - >output.sam &&
        cat output.sam | awk '{print $3}' | sort -k1,1V | uniq -c >count_no_mapq0.txt &&
        awk '! /random/ && ! /Un/ && /chr/  ' count_no_mapq0.txt | awk '{print $2, $1}' OFS="\t" | sort -k1,1 -V -s >temp2.2.txt

    if [[ $types == SE ]]; then
        mv temp2.2.txt 'step2.2_chrom_count_'$name'.txt'
    elif [[ $types == PE ]]; then
        awk '{$2=int($2*0.5); print}' OFS="\t" temp2.2.txt >'step2.2_chrom_count_'$name'.txt' && rm temp2.2.txt
    fi

    # only effect reads
    $methylQA density -S $chrom_size output.sam
    cut -f 1 output.extended.bed | uniq -c >count_unique.txt
    awk '{print $2, $1}' OFS="\t" count_unique.txt | sort -k1,1 -V -s >'step2.2_chrom_count_unique_'$name'.txt'
    effect_chrM=$(grep chrM output.extended.bed | wc -l)

    # get chrM count in uniquely mapped reads
    # to keep the count consistent, the results are all from methylQA
    # when count directly from output.sam file, please pay attention to the unpaired reads
    cat <(samtools view -H input.sam) <(awk '$3=="chrM"' input.sam) | $methylQA density -S -r -o temp $chrom_size -
    unique_chrM=$(grep chrM temp.extended.bed | wc -l)

    rm temp*
    rm output*
    rm count*.txt
    rm input.sam
    awk -F "\t" '{print $2}' 'step2.2_chrom_count_unique_'$name'.txt' | paste 'step2.2_chrom_count_'$name'.txt' - | awk -F "\t" -v marker=$marker '{print $1,$2+0,$3+0,marker}' OFS="\t" | awk 'length($1)<7' >./'QC_data_collection_'$name/'step2.2_chrom_count_'$name'.result'

    if [ $? == 0 ]; then
        echo "step2.2, count reads distribution process done" >>QC_pipe_processing.log
    else
        echo "step2.2, count reads distribution process fail......" >>QC_pipe_processing.log
    fi

    rm step2.2_chrom_count*txt
}
```

**2.3 Library complexity estimate**

Tool: preseq v2.0.0
Input: aligned BAM file
```
$preseq lc_extrap -o 'step2.3_yield_'$name'.result' -B 'step2.1_trimed_'$name'.bam'
```

### Step3, BAM file processing

Tool: methylQA v0.2.1, MACS2 2.2.7.1<br/>
Input: aligned BAM file, reads BED file
```
s3.1_methylQA() {
    echo 'step3.1 methylQA processing......'
    grep -v ^chrM $chrom_size >nochrM_chrom_size.txt &&
        rm $chrom_size &&
        chrom_size=nochrM_chrom_size.txt
    echo "step3.1, the chrom file for methylQA is $chrom_size" >>QC_pipe_processing.log

    if [[ $data == ATAC ]]; then
        echo "methylQA min insertion length choice $methylQA_cutoff" >>QC_pipe_processing.log
        $methylQA atac -X $methylQA_cutoff -o step3.1_methylQA_$name $chrom_size 'step2.1_trimed_'$name'.bam'
        if [ $? == 0 ]; then
            echo "step3.1, mathylQA atac process done" >>QC_pipe_processing.log
        else
            echo "step3.1, mathylQA atac process fail......" >>QC_pipe_processing.log
        fi
    elif [[ $data == CUTnTag ]]; then
        echo "step3.1, there is no need for CUT&Tag-seq data having the methylQA min insertion length" >>QC_pipe_processing.log
        $methylQA density -o step3.1_methylQA_$name $chrom_size 'step2.1_trimed_'$name'.bam'
        $methylQA density -o step3.1_methylQA_$control_name $chrom_size 'step2.1_trimed_'$control_name'.bam'
        if [ $? == 0 ]; then
            echo "step3.1, mathylQA density process done" >>QC_pipe_processing.log
        else
            echo "step3.1, mathylQA density process fail......" >>QC_pipe_processing.log
        fi
    fi

    # mapping status
    map_mapped=$(grep 'mappable reads' step3.1_methylQA_$name'.report' | awk '{print $4}')
    map_uniq=$(grep '(mapQ >= 10)' step3.1_methylQA_$name'.report' | awk '{print $8}')
    map_effect=$(grep 'non-redundant' step3.1_methylQA_$name'.report' | awk '{print $6}')
    mapped_ratio=$(echo "scale=2; ($map_mapped/$raw_reads)" | bc -l)
    effect_ratio=$(echo "scale=2; ($map_effect)/($raw_reads)" | bc -l)

    # unique chrM ratio from step2.2
    unique_chrM_ratio=$(python3 -c "print($unique_chrM*1.0 / ($unique_chrM+$map_uniq) )")
    echo -e "non_chrM_unique_mapped\tchrM\tunique_chrM_ratio" >'step2.2_unique_chrM_ratio_'$name'.result'
    echo -e "$map_uniq\t$unique_chrM\t$unique_chrM_ratio" >>'step2.2_unique_chrM_ratio_'$name'.result'
    mv 'step2.2_unique_chrM_ratio_'$name'.result' ./'QC_data_collection_'"$name"

    #unique_no_chrM=`python3 -c "print($map_uniq-$unique_chrM)"`
    #effect_no_chrM=`python3 -c "print($map_effect-$effect_chrM)"`
    nodup_ratio=$(echo "scale=3; $map_effect/$map_uniq" | bc -l)
    after_dup=$(python3 -c "print(1-$nodup_ratio*1.0)")

    useful=$(grep 'non-redundant' step3.1_methylQA_"$name"*.report | awk '{print $6}')
    if [[ $data == ATAC ]]; then
        single_end=$(wc -l *open.bed | awk '{print $1}')
    elif [[ $data == CUTnTag ]]; then
        single_end=$(wc -l *"$name"*extended.bed | awk '{print $1}')
    fi
    uf_ratio=$(echo "scale=3; $useful / $raw_reads" | bc -l)
    echo -e "file\ttotal\tuseful\tuseful_ratio\tsingle_end" >'step3.1_useful_reads_'$name.result
    echo -e "$name\t$raw_reads\t$useful\t$uf_ratio\t$single_end" >>'step3.1_useful_reads_'$name.result
    mv 'step3.1_useful_reads_'$name.result ./'QC_data_collection_'$name
    if [[ $data == ATAC ]]; then
        sort -n 'step3.1_methylQA_'*$name'.insertdistro' | uniq -c | awk '{print $2,$1}' >'step3.1_insertion_distri_'$name'.result' && rm 'step3.1_methylQA_'*$name'.insertdistro'
        mv 'step3.1_insertion_distri_'$name'.result' ./'QC_data_collection_'$name
        rm step3.1_methylQA_*bigWig
    elif [[ $data == CUTnTag ]]; then
        echo -e "step3.1_distri: there is no insertion distribution plot for CUT&Tag-seq data, skip it" >>QC_pipe_processing.log
        mv 'step3.1_methylQA_'$name'.genomeCoverage' ./'QC_data_collection_'$name
        mv 'step3.1_methylQA_'$control_name'.genomeCoverage' ./'QC_data_collection_'$name
    fi
    mv *pdf ./'QC_data_collection_'$name
}

# 3.2, blacklist-removal and normalization
s3.2_rmbl_bg() {
    echo 'step3.2 blacklist-removal and normalization'

    if [[ $data == ATAC ]]; then
        echo 'remove black list, and do normalization(only for ATAC-seq data)...'
        # add a new bigwig file without black list
        intersectBed -iobuf 200M -a 'step3.1_methylQA_'*$name'.open.bedGraph' -b $black_list -v >rmbl.bedGraph
        bedGraphToBigWig rmbl.bedGraph $chrom_size 'step3.2_rmbl_'$name'.bigWig' && rm 'step3.1_methylQA_'*$name'.open.bedGraph'

        # normalization
        norm=$(grep 'non-redundant' step3.1_methylQA_*report | awk '{print $6}')
        factor=$(echo "scale=3; $norm/10000000" | bc -l)
        awk -v factor=$factor '{print $1,$2,$3,$4/factor}' OFS='\t' rmbl.bedGraph >'step3.2_normalized_per_10M_'$name'.open.bedGraph'
        bedGraphToBigWig 'step3.2_normalized_per_10M_'$name'.open.bedGraph' $chrom_size 'step3.2_normalized_per_10M_'$name'.bigWig' && rm 'step3.2_normalized_per_10M_'$name'.open.bedGraph'
        rm rmbl.bedGraph # because for ATAC-seq data, we could not use SEACR to call peaks, so there is no need to keep bedGrapg file

        if [ $? == 0 ]; then
            echo "step3.2, ATAC-seq normalization process done" >>QC_pipe_processing.log
        else
            echo "step3.2, ATAC-seq normalization process failed......" >>QC_pipe_processing.log
        fi

    elif [[ $data == CUTnTag ]]; then
        echo "step3.2, CUT&Tag: we remove black list to create new bedGraph and BigWig files, but normalization process is not necessary for CUT&Tag-seq data, just skip" >>QC_pipe_processing.log
        # for Marker_input data
        intersectBed -iobuf 200M -a 'step3.1_methylQA_'$name'.extended.bedGraph' -b $black_list -v >'step3.2_rmbl_'"$name"'.bedGraph'
        bedGraphToBigWig 'step3.2_rmbl_'"$name"'.bedGraph' $chrom_size 'step3.2_rmbl_'$name'.bigWig' && rm 'step3.1_methylQA_'*$name'.extended.bedGraph'
        # for IgG_control data
        intersectBed -iobuf 200M -a 'step3.1_methylQA_'$control_name'.extended.bedGraph' -b $black_list -v >'step3.2_rmbl_'"$control_name"'.bedGraph'
        bedGraphToBigWig 'step3.2_rmbl_'"$control_name"'.bedGraph' $chrom_size 'step3.2_rmbl_'$control_name'.bigWig' && rm 'step3.1_methylQA_'*$control_name'.extended.bedGraph'
        mv 'step3.1_methylQA_'$control_name'.bigWig' ./'QC_data_collection_'$name
        if [ $? == 0 ]; then
            echo "step3.2, CUT&Tag-seq blacklist-removal process done" >>QC_pipe_processing.log
        else
            echo "step3.2, CUT&Tag-seq blacklist-removal process failed......" >>QC_pipe_processing.log
        fi
    fi

    mv *report ./'QC_data_collection_'$name 2>/dev/null
}
```
note: methylQA atac would do a PE->SE transformation based on insertion site for PE data.


### Step4, Peak calling

**4.0 IgG data preparation**
```
s4.0_IgG_data_pre() {
    echo 'step4.0 IgG_data preparation'

    if [[ $data == ATAC ]]; then
        echo 'step4.0, there is no IgG data for ATAC-seq data, skip this step......' >>QC_pipe_processing.log

    elif [[ $data == CUTnTag ]]; then
        mkdir 'QC_IgG_data_collection_'"$control_name"
        mv ./'QC_data_collection_'"$name"/*"$control_name"* ./'QC_IgG_data_collection_'"$control_name"
        mv ./'step'*"$control_name"* ./'QC_IgG_data_collection_'"$control_name"

        if [ $? == 0 ]; then
            echo "step4.0, CUT&Tag-seq IgG data prepare process done" >>QC_pipe_processing.log
        else
            echo "step4.0, CUT&Tag-seq IgG data prepare process failed......" >>QC_pipe_processing.log
        fi
    fi
}
```

**4.1 Peak calling**
Tool: Peak calling software defined in command line
Input: software required input files
Output: BED format file containing peaks
```
s4.1_peakcall() {
    echo 'step4.1 peakcalling...'

    if [[ $data == ATAC ]]; then
        echo 'peak calling step for ATAC-seq data...'
        if [[ $software == AIAP ]]; then
            echo 'use AIAP (optimal MACS2 parameters) to call peak'
            # prepare input BED file for MACS2
            awk '{if ($2 > $3)sub($2, 0); print}' OFS="\t" 'step3.1_methylQA_'$name'.open.bed' >temp.open.bed &&
                intersectBed -iobuf 200M -a temp.open.bed -b $black_list -v >'step3.3_rmbl_'$name'.open.bed' &&
                rm temp.open.bed 'step3.1_methylQA_'$name'.open.bed'
            # MACS2 call peak
            if [[ $personalize == False ]]; then
                $macs2 callpeak -t 'step3.3_rmbl_'$name'.open.bed' -g "$macs2_genome" -q 0.01 -n 'step4.1_peakcall_AIAP_'$name --keep-dup 1000 --nomodel --shift 0 --extsize 150
                # peak length distribution:
                awk '{print $3-$2+1}' 'step4.1_peakcall_AIAP_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                mv 'step4.1_peakcall_AIAP_'$name'_peaks.xls' ./'QC_data_collection_'$name
                mv 'step4.1_peakcall_AIAP_'$name'_summits.bed' ./'QC_data_collection_'$name
            elif [[ $personalize == True ]]; then
                if [[ $callsummits == False ]]; then
                    if [[ $nomodel == False ]] && [[ $broad == False ]]; then
                        $macs2 callpeak -t 'step3.3_rmbl_'$name'.open.bed' -g "$macs2_genome" -n 'step4.1_peakcall_AIAP_'$name -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_AIAP_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_summits.bed' ./'QC_data_collection_'$name
                        mv ./*'model.r' ./'QC_data_collection_'$name
                    elif [[ $nomodel == True ]] && [[ $broad == False ]]; then
                        $macs2 callpeak -t 'step3.3_rmbl_'$name'.open.bed' -g "$macs2_genome" -n 'step4.1_peakcall_AIAP_'$name -q $qvalue --keep-dup $keepdup --nomodel --shift $shift --extsize $extsize
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_AIAP_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_summits.bed' ./'QC_data_collection_'$name
                    elif [[ $nomodel == True ]] && [[ $broad == True ]]; then
                        $macs2 callpeak -t 'step3.3_rmbl_'$name'.open.bed' -g "$macs2_genome" -n 'step4.1_peakcall_AIAP_'$name -q $qvalue --keep-dup $keepdup --nomodel --shift $shift --extsize $extsize --broad
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_AIAP_'$name'_peaks.broadPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_peaks.gappedPeak' ./'QC_data_collection_'$name
                    elif [[ $nomodel == False ]] && [[ $broad == True ]]; then
                        $macs2 callpeak -t 'step3.3_rmbl_'$name'.open.bed' -g "$macs2_genome" -n 'step4.1_peakcall_AIAP_'$name -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --broad
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_AIAP_'$name'_peaks.broadPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_peaks.gappedPeak' ./'QC_data_collection_'$name
                        mv ./*'model.r' ./'QC_data_collection_'$name
                    fi
                elif [[ $callsummits == True ]]; then
                    if [[ $nomodel == False ]] && [[ $broad == False ]]; then
                        $macs2 callpeak -t 'step3.3_rmbl_'$name'.open.bed' -g "$macs2_genome" -n 'step4.1_peakcall_AIAP_'$name -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --call-summits
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_AIAP_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_summits.bed' ./'QC_data_collection_'$name
                        mv ./*'model.r' ./'QC_data_collection_'$name
                    elif [[ $nomodel == True ]] && [[ $broad == False ]]; then
                        $macs2 callpeak -t 'step3.3_rmbl_'$name'.open.bed' -g "$macs2_genome" -n 'step4.1_peakcall_AIAP_'$name -q $qvalue --keep-dup $keepdup --nomodel --shift $shift --extsize $extsize --call-summits
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_AIAP_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_AIAP_'$name'_summits.bed' ./'QC_data_collection_'$name
                    elif [[ $nomodel == True ]] && [[ $broad == True ]]; then
                        echo "--call-summits and --broad could NOT be used at the same time..."
                        exit
                    elif [[ $nomodel == False ]] && [[ $broad == True ]]; then
                        echo "--call-summits and --broad could NOT be used at the same time..."
                        exit
                    fi
                fi
            fi

        elif [[ $software == F-seq ]]; then
            # prepare input BED file for F-seq
            awk '{if ($2 > $3)sub($2, 0); print}' OFS="\t" 'step3.1_methylQA_'$name'.open.bed' >temp.open.bed &&
                intersectBed -iobuf 200M -a temp.open.bed -b $black_list -v >'step3.3_rmbl_'$name'.open.bed' &&
                rm temp.open.bed 'step3.1_methylQA_'$name'.open.bed'
            # F-seq call peak
            $fseq -of bed 'step3.3_rmbl_'$name'.open.bed'
            cat chr* >'step4.1_peakcall_fseq_'"$name"'.bed'
            rm chr*
            # peak length distribution:
            awk '{print $3-$2+1}' 'step4.1_peakcall_fseq_'"$name"'.bed' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
            mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
            # mv step3.1_methylQA_"$name"*  ./'QC_data_collection_'"$name"
        elif [[ $software == HMMRATAC ]]; then
            if [[ $types == PE ]]; then
                # prepare input BED file for MACS2
                awk '{if ($2 > $3)sub($2, 0); print}' OFS="\t" 'step3.1_methylQA_'$name'.open.bed' >temp.open.bed &&
                intersectBed -iobuf 200M -a temp.open.bed -b $black_list -v >'step3.3_rmbl_'$name'.open.bed' &&
                rm temp.open.bed 'step3.1_methylQA_'$name'.open.bed'
                #prepare the required files
                $samtools sort 'step2.1_trimed_'"$name"'.bam' >'temp_'"$name"'.sorted.bam'
                $samtools index 'temp_'"$name"'.sorted.bam' >'temp_'"$name"'.sorted.bam.bai'
                $samtools view -H 'temp_'"$name"'.sorted.bam' | $perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' >'temp_'"$name"'.genome.info'
                if [[ $personalize == False ]]; then
                    # call peak by HMMRATAC
                    $java -jar $hmmr_path -b 'temp_'"$name"'.sorted.bam' -i 'temp_'"$name"'.sorted.bam.bai' -g 'temp_'"$name"'.genome.info' -o 'step4.1_peakcall_HMMR_'$name --bedgraph True -u 35 -l 2
                elif [[ $personalize == True ]]; then
                    $java -jar $hmmr_path -b 'temp_'"$name"'.sorted.bam' -i 'temp_'"$name"'.sorted.bam.bai' -g 'temp_'"$name"'.genome.info' -o 'step4.1_peakcall_HMMR_'$name --bedgraph $bedgraph -u $upper -l $lower -e $blacklist --score $score
                fi
                # create BED3 format output
                cut -f 1-3 'step4.1_peakcall_HMMR_'$name'_peaks.gappedPeak' | sort -k1,1 -k2,2n >'step4.1_peakcall_HMMR_BED3_'$name'.bed'
                # peak length distribution:
                awk '{print $3-$2+1}' 'step4.1_peakcall_HMMR_BED3_'$name'.bed' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                # remove temp files
                rm temp*
                mv 'step4.1_peakcall_HMMR_'$name* ./'QC_data_collection_'$name

            elif [[ $types == SE ]]; then
                echo 'Error: step4.1, HMMRATAC is designed for Pair-end ATAC-seq data, it could not run for Single-end data' >>QC_pipe_processing.log
            fi

        else
            echo 'step4.1, the software you select is not suitable for ATAC-seq data peakcalling, please select again' >>QC_pipe_processing.log
        fi

    elif [[ $data == CUTnTag ]]; then
        echo 'peak calling step for CUT&Tag-seq data...'

        if [[ $software == MACS2 ]]; then
            if [[ $marker == narrow ]]; then
                echo 'use MACS2 BED IgG to do narrow peakcalling'
                # MACS2 call peak
                if [[ $personalize == False ]]; then
                    $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -g $macs2_genome -q 0.01 -n 'step4.1_peakcall_MACS2_'$name
                    # peak length distribution:
                    awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                    mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                    mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                    mv 'step4.1_peakcall_MACS2_'$name'_summits.bed' ./'QC_data_collection_'$name
                    mv ./*'model.r' ./'QC_data_collection_'$name
                elif [[ $personalize == True ]]; then
                    if [[ $callsummits == False ]]; then
                        if [[ $nomodel == False ]] && [[ $broad == False ]]; then
                            $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize
                            # peak length distribution:
                            awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                            mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_summits.bed' ./'QC_data_collection_'$name
                            mv ./*'model.r' ./'QC_data_collection_'$name 2>/dev/null
                        elif [[ $nomodel == False ]] && [[ $broad == True ]]; then
                            $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --broad
                            # peak length distribution:
                            awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.broadPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                            mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_peaks.gappedPeak' ./'QC_data_collection_'$name
                            mv ./*'model.r' ./'QC_data_collection_'$name 2>/dev/null
                        elif [[ $nomodel == True ]] && [[ $broad == True ]]; then
                            $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --nomodel --broad
                            # peak length distribution:
                            awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.broadPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                            mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_peaks.gappedPeak' ./'QC_data_collection_'$name
                        elif [[ $nomodel == True ]] && [[ $broad == False ]]; then
                            $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --nomodel
                            # peak length distribution:
                            awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                            mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_summits.bed' ./'QC_data_collection_'$name
                        fi
                    elif [[ $callsummits == True ]]; then
                        if [[ $nomodel == False ]] && [[ $broad == False ]]; then
                            $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --call-summits
                            # peak length distribution:
                            awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                            mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_summits.bed' ./'QC_data_collection_'$name
                            mv ./*'model.r' ./'QC_data_collection_'$name 2>/dev/null
                        elif [[ $nomodel == False ]] && [[ $broad == True ]]; then
                            echo "--call-summits and --broad could NOT be used at the same time..."
                            exit
                        elif [[ $nomodel == True ]] && [[ $broad == True ]]; then
                            echo "--call-summits and --broad could NOT be used at the same time..."
                            exit
                        elif [[ $nomodel == True ]] && [[ $broad == False ]]; then
                            $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --nomodel --call-summits
                            # peak length distribution:
                            awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                            mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                            mv 'step4.1_peakcall_MACS2_'$name'_summits.bed' ./'QC_data_collection_'$name
                        fi
                    fi
                fi
            elif [[ $marker == broad ]]; then
                echo 'use MACS2 BED IgG to do broad peakcalling'
                # MACS2 call peak
                if [[ $personalize == False ]]; then
                    $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -g $macs2_genome -q 0.01 -n 'step4.1_peakcall_MACS2_'$name --broad
                    # peak length distribution:
                    awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.broadPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                    mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                    mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                    mv 'step4.1_peakcall_MACS2_'$name'_peaks.gappedPeak' ./'QC_data_collection_'$name
                    mv ./*'model.r' ./'QC_data_collection_'$name 2>/dev/null
                elif [[ $personalize == True ]]; then
                    if [[ $nomodel == False ]] && [[ $broad == False ]]; then
                        $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_MACS2_'$name'_summits.bed' ./'QC_data_collection_'$name
                        mv ./*'model.r' ./'QC_data_collection_'$name 2>/dev/null
                    elif [[ $nomodel == False ]] && [[ $broad == True ]]; then
                        $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --broad
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.broadPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_MACS2_'$name'_peaks.gappedPeak' ./'QC_data_collection_'$name
                        mv ./*'model.r' ./'QC_data_collection_'$name 2>/dev/null
                    elif [[ $nomodel == True ]] && [[ $broad == True ]]; then
                        $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --nomodel --broad
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.broadPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_MACS2_'$name'_peaks.gappedPeak' ./'QC_data_collection_'$name
                    elif [[ $nomodel == True ]] && [[ $broad == False ]]; then
                        $macs2 callpeak -t 'step3.1_methylQA_'$name'.extended.bed' -c ./'QC_IgG_data_collection_'"$control_name"/'step3.1_methylQA_'"$control_name"'.extended.bed' -n 'step4.1_peakcall_MACS2_'$name -g $macs2_genome -q $qvalue --keep-dup $keepdup --shift $shift --extsize $extsize --nomodel
                        # peak length distribution:
                        awk '{print $3-$2+1}' 'step4.1_peakcall_MACS2_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                        mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_MACS2_'$name'_peaks.xls' ./'QC_data_collection_'$name
                        mv 'step4.1_peakcall_MACS2_'$name'_summits.bed' ./'QC_data_collection_'$name
                    fi
                fi
            fi
        elif [[ $software == F-seq ]]; then
            # F-seq call peak
            $fseq -of bed 'step3.1_methylQA_'$name'.extended.bed'
            cat chr* >'step4.1_peakcall_fseq_'"$name"'.bed'
            rm chr*
            # peak length distribution:
            awk '{print $3-$2+1}' 'step4.1_peakcall_fseq_'"$name"'.bed' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
            mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name

        elif [[ $software == SEACR ]]; then
            # SEACR call peak
            if [[ $personalize == False ]]; then
                bash $seacr_path 'step3.2_rmbl_'$name'.bedGraph' ./'QC_IgG_data_collection_'"$control_name"/'step3.2_rmbl_'"$control_name"'.bedGraph' norm stringent $name
                # peak length distribution:
                awk '{print $3-$2+1}' $name'.stringent.bed' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                mv $name'.stringent.bed' 'step4.1_peakcall_SEACR_'"$name"'.bed'
            elif [[ $personalize == True ]]; then
                if [[ $model == stringent ]]; then
                    bash $seacr_path 'step3.2_rmbl_'$name'.bedGraph' ./'QC_IgG_data_collection_'"$control_name"/'step3.2_rmbl_'"$control_name"'.bedGraph' $normalize stringent $name
                    # peak length distribution:
                    awk '{print $3-$2+1}' $name'.stringent.bed' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                    mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                    mv $name'.stringent.bed' 'step4.1_peakcall_SEACR_'"$name"'.bed'
                elif [[ $model == relaxed ]]; then
                    bash $seacr_path 'step3.2_rmbl_'$name'.bedGraph' ./'QC_IgG_data_collection_'"$control_name"/'step3.2_rmbl_'"$control_name"'.bedGraph' $normalize relaxed $name
                    # peak length distribution:
                    awk '{print $3-$2+1}' $name'.relaxed.bed' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                    mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                    mv $name'.relaxed.bed' 'step4.1_peakcall_SEACR_'"$name"'.bed'
                fi
            fi

        elif [[ $software == HMMRATAC ]]; then
            if [[ $types == PE ]]; then
                #prepare the required files
                $samtools sort 'step2.1_trimed_'"$name"'.bam' >'temp_'"$name"'.sorted.bam'
                $samtools index 'temp_'"$name"'.sorted.bam' >'temp_'"$name"'.sorted.bam.bai'
                $samtools view -H 'temp_'"$name"'.sorted.bam' | $perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' >'temp_'"$name"'.genome.info'
                # call peak by HMMRATAC
                $java -jar $hmmr_path -b 'temp_'"$name"'.sorted.bam' -i 'temp_'"$name"'.sorted.bam.bai' -g 'temp_'"$name"'.genome.info' -o 'step4.1_peakcall_HMMR_'$name --bedgraph True -u 35 -l 2
                # create BED3 format output
                cut -f 1-3 'step4.1_peakcall_HMMR_'$name'_peaks.gappedPeak' | sort -k1,1 -k2,2n >'step4.1_peakcall_HMMR_BED3_'$name'.bed'
                # peak length distribution:
                awk '{print $3-$2+1}' 'step4.1_peakcall_HMMR_BED3_'$name'.bed' | sort -n | uniq -c | awk '{print $2,$1}' >'step4.1_peak_length_distri_'$name'.result'
                mv 'step4.1_peak_length_distri_'$name'.result' ./'QC_data_collection_'$name
                # remove temp files
                rm temp*
                mv 'step4.1_peakcall_HMMR_'$name* ./'QC_data_collection_'$name
            elif [[ $types == SE ]]; then
                echo 'Error: step4.1, HMMRATAC is designed for Pair-end data, it could not run for Single-end data' >>QC_pipe_processing.log
            fi

        elif [[ $software == ChromHMM ]]; then
            # prepare needed directories and files
            mkdir ./'step4.1_peakcall_ChromHMM_'"$name"'_dir'
            mkdir ./'step4.1_peakcall_ChromHMM_'"$name"'_dir'/'temp_binary_output_dir'
            mkdir ./'step4.1_peakcall_ChromHMM_'"$name"'_dir'/'input_data_dir'
            mkdir ./'step4.1_peakcall_ChromHMM_'"$name"'_dir'/'output_dir'
            cp 'step2.1_trimed_'"$name"'.bam' ./'step4.1_peakcall_ChromHMM_'"$name"'_dir'/'input_data_dir'
            cp ./'QC_IgG_data_collection_'"$control_name"/'step2.1_trimed_'"$control_name"'.bam' ./'step4.1_peakcall_ChromHMM_'"$name"'_dir'/'input_data_dir'
            cd ./'step4.1_peakcall_ChromHMM_'"$name"'_dir'/'input_data_dir' || exit
            Marker_bam='step2.1_trimed_'"$name"'.bam'
            IgG_bam='step2.1_trimed_'"$control_name"'.bam'
            echo "$species" "$name" "$Marker_bam" "$IgG_bam" >'temp_input_table_'"$name"'.txt'
            awk -v FS=' ' -v OFS='\t' '{print $1,$2,$3,$4}' 'temp_input_table_'"$name"'.txt' >'input_table_'"$name"'.txt'
            rm 'temp_input_table_'"$name"'.txt'
            cd ..
            # call peak by ChromHMM
            if [[ $personalize == False ]]; then
                # use BinarizeBam to train data (generate binarized files)
                $java -mx4000M -jar "$chromhmm_path" BinarizeBam ../'nochrM_chrom_size.txt' ./'input_data_dir' ./'input_data_dir'/'input_table_'"$name"'.txt' ./'temp_binary_output_dir'
                # use LearnModel to get segment results
                $java -mx4000M -jar "$chromhmm_path" LearnModel ./'temp_binary_output_dir' ./'output_dir' 2 "$species"
                # select annotated regions
                awk '{if($4==2) print $0}' ./'output_dir'/"$species"'_2_dense.bed' | cut -f 1-3 | sort -k1,1 -k2,2n >./'output_dir'/'step4.1_peakcall_ChromHMM_BED3_'"$name"'.bed'
                cp ./'output_dir'/'step4.1_peakcall_ChromHMM_BED3_'"$name"'.bed' ../
                rm -rf ./'temp_binary_output_dir'
                cd ..
                mv ./'step4.1_peakcall_ChromHMM_'"$name"'_dir' ./'QC_data_collection_'"$name"
            elif [[ $personalize == True ]]; then
                # use BinarizeBam to train data (generate binarized files)
                $java -mx4000M -jar "$chromhmm_path" BinarizeBam -b $binsize -f $foldthresh ../'nochrM_chrom_size.txt' ./'input_data_dir' ./'input_data_dir'/'input_table_'"$name"'.txt' ./'temp_binary_output_dir'
                # use LearnModel to get segment results
                $java -mx4000M -jar "$chromhmm_path" LearnModel -b $binsize -init $init -s $seed -z $zerotransitionpower ./'temp_binary_output_dir' ./'output_dir' 2 "$species"
                # select annotated regions
                awk '{if($4==2) print $0}' ./'output_dir'/*'_dense.bed' | cut -f 1-3 | sort -k1,1 -k2,2n >./'output_dir'/'step4.1_peakcall_ChromHMM_BED3_'"$name"'.bed'
                cp ./'output_dir'/'step4.1_peakcall_ChromHMM_BED3_'"$name"'.bed' ../
                rm -rf ./'temp_binary_output_dir'
                cd ..
                mv ./'step4.1_peakcall_ChromHMM_'"$name"'_dir' ./'QC_data_collection_'"$name"
            fi
        fi
    fi

    if [ $? == 0 ]; then
        echo "step4.1, peakcalling process done" >>QC_pipe_processing.log
    else
        echo "step4.1, peakcalling prepare process failed......" >>QC_pipe_processing.log
    fi

}
```




































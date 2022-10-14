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












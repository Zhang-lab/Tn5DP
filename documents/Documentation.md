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














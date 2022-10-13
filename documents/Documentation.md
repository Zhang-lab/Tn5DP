# Documentation v1.0
ATAC-seq and CUT&Tag-seq data analysis pipeline for Bo Zhang`s lab<br/>Last edit: 10/12/2022<br/>For any question please contact: siyuancheng@wustl.edu

**Outline**
1. Parameters and options
2. Output example and annotation
3. Data processing details

<br />

## 1. Parameters and options
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













# Tn5DP (Tn5-based epigenetic Data Processing)
Pipeline for the QC metrics construction, data analysis and visualization of ATAC-seq and CUT&Tag-seq data. 

Current Version: `Tn5DP_v1.0` Last update: `2022.10.12`


Advisor: Bo Zhang<br/>Contributor: Siyuan Cheng

For any question, please contact [Siyuan Cheng](siyuancheng@wustl.ed):point_left:

<br />
<br /> 

## Documentation
1. Pipeline documentation: analysis details and QC metrics information<br/>Please [click here](documents/Documentation.md)
2. Update logfile: pipeline change record<br/>Please [click here](documents/update_log.md)

<br />
<br />

## Usage:
### Singularity3 Installation
Tn5DP image has to be run with **Singularity version3+**, you could follow this instruction if you haven`t install Singularity3. <br/>Please [click here](https://github.com/sylabs/singularity/blob/main/INSTALL.md)<br/>(You will need sudo permission to properlly install and configure it, but you can run it without sudo after installation:smiley:)

### Test ATAC-seq data:
There are one paired-end mm10 data with 0.25M reads for test purpose, they can be downloaded by:
```
wget http://regmedsrv1.wustl.edu/Public_SPACE/ryan/Public_html/singularity_ac/sample_data/atac-seq/mm10_1.fastq.gz
wget http://regmedsrv1.wustl.edu/Public_SPACE/ryan/Public_html/singularity_ac/sample_data/atac-seq/mm10_2.fastq.gz
```

### Test CUT&Tag-seq data:
There are one paired-end mm10 data (H3K36me3 marker) with 0.25M reads for test purpose, they can be downloaded by:
```
wget http://regmedsrv1.wustl.edu/Public_SPACE/ryan/Public_html/singularity_ac/sample_data/cuttag-seq/mm10_k36me3_1.fastq.gz
wget http://regmedsrv1.wustl.edu/Public_SPACE/ryan/Public_html/singularity_ac/sample_data/cuttag-seq/mm10_k36me3_2.fastq.gz
wget http://regmedsrv1.wustl.edu/Public_SPACE/ryan/Public_html/singularity_ac/sample_data/cuttag-seq/mm10_igg_1.fastq.gz
wget http://regmedsrv1.wustl.edu/Public_SPACE/ryan/Public_html/singularity_ac/sample_data/cuttag-seq/mm10_igg_2.fastq.gz
```












#!/bin/bash
# This is for Docker pipe, default root dir is /atac-seq
# This is for Docker!!!

species=$1

# get pipe path, though readlink/realpath can do it, some version doesn't have that
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# common tools:
adapter_1="CTGTCTCTTATACACATCT"
adapter_2="CTGTCTCTTATACACATCT"
fastq_dump_tool='/usr/local/ncbi/sra-tools/bin/fastq-dump'
preseq="/usr/bin/preseq"
cutadapt="/home/ryan/.local/bin/cutadapt"
fastqc="/usr/bin/fastqc"
samtools="/usr/local/bin/samtools"
bwa="/usr/bin/bwa"
methylQA="/usr/bin/methylQA"
macs2="/home/ryan/Python-3.9.6/bin/macs2"
# added tools which are not included by AIAP
java="/home/ryan/java1.8/jdk1.8/bin/java"
perl="/usr/bin/perl"
fseq="/home/ryan/controlsoftware/F-seq/dist~/fseq/bin/fseq"
hmmr_path="/home/ryan/controlsoftware/hmmr/HMMRATAC-master/HMMRATAC_V1.2.10_exe.jar"
seacr_path="/home/ryan/controlsoftware/SEACR-master/SEACR_1.3.sh"
chromhmm_path="/home/ryan/controlsoftware/ChromHMM/ChromHMM.jar"



# genome specific resources:
if [[ $species == mm39 ]]; 
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/mm39/mm39_bwa.index/mm39.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/mm39/mm39.chrom.sizes"
	black_list="/home/ryan/benchmark/singularity/Genome/mm39/mm39_black_list.bed"
	genome_size=2728222451
	promoter_file="/home/ryan/benchmark/singularity/Genome/mm39/mm39_1kb_promoter.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/mm39/mm39_1kb_coding_promoter.bed"
	macs2_genome='mm'
elif [[ $species == mm10 ]]; 
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/mm10/bwa_index_mm10/mm10.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/mm10/mm10.chrom.sizes"
	black_list="/home/ryan/benchmark/singularity/Genome/mm10/mm10_black_list.bed"
	genome_size=2730871774
	promoter_file="/home/ryan/benchmark/singularity/Genome/mm10/mm10_promoter_bistream_1kb.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/mm10/mm10_promoter_coding_bistream_1kb.bed"
	macs2_genome='mm'
elif [[ $species == mm9 ]];
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/mm9/bwa_index_mm9/mm9.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/mm9/mm9_chrom_sizes"
	black_list="/home/ryan/benchmark/singularity/Genome/mm9/mm9-blacklist.bed"
	genome_size=2725765481
	promoter_file="/home/ryan/benchmark/singularity/Genome/mm9/mm9_promoter_bistream_1kb.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/mm9/mm9_coding_promoter_bistream_1kb.bed"
	macs2_genome='mm'
elif [[ $species == hg38 ]];
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/hg38/bwa_index_hg38.25/hg38.25_chromsome.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/hg38/hg38.25_chromsome.sizes"
	black_list="/home/ryan/benchmark/singularity/Genome/hg38/hg38_black_list.bed"
	genome_size=3209286105
	promoter_file="/home/ryan/benchmark/singularity/Genome/hg38/hg38_promoter_bistream_1kb.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/hg38/hg38_coding_promoter_bistream_1kb.bed"
	macs2_genome='hs'
elif [[ $species == hg19 ]];
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/hg19/bwa_index_0.7.5/hg19.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/hg19/hg19_chromosome.size"
	black_list="/home/ryan/benchmark/singularity/Genome/hg19/hg19_blacklist.bed"
	genome_size=3137161264
	promoter_file="/home/ryan/benchmark/singularity/Genome/hg19/hg19_promoter_bistream_1kb.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/hg19/hg19_promoter_coding_bistream_1kb.bed"
	macs2_genome='hs'
elif [[ $species == danRer10 ]];
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/danRer10/bwa_index_denRer10/danRer10.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/danRer10/danRer10.chrom.sizes"
	touch pesudo_bl.txt
	black_list="pesudo_bl.txt"
	genome_size=1340447187
	promoter_file="/home/ryan/benchmark/singularity/Genome/danRer10/promoter_region_danRer10_bistream_1k.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/danRer10/danRer10_coding_promoter_bistream_1k.bed"
	macs2_genome='mm'
elif [[ $species == danRer11 ]];
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/danRer11/BWAIndex/danRer11.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/danRer11/GRCz11_chrom.size"
	touch pesudo_bl.txt
    black_list="pesudo_bl.txt"
    genome_size=1345118429
	promoter_file="/home/ryan/benchmark/singularity/Genome/danRer11/GRCz11_promoter_region.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/danRer11/pseudo_GRCz11_coding_promoter_region.bed"
	macs2_genome='mm'
elif [[ $species == dm6 ]];
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/dm6/bwa_index_dm6/d.mel.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/dm6/d.mel.chrom.sizes"
	touch pesudo_bl.txt
    black_list="pesudo_bl.txt"
    genome_size=143726002
	promoter_file="/home/ryan/benchmark/singularity/Genome/dm6/promoter_region_from_Dmel.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/dm6/pseudo_coding_promoter_region.bed"
	macs2_genome="dm"
elif [[ $species == rn6 ]];
	then
	bwa_ref="/home/ryan/benchmark/singularity/Genome/rn6/bwa_index_rn6/rn6.fa"
	chrom_size="/home/ryan/benchmark/singularity/Genome/rn6/rn6.chrom.sizes"
	touch pesudo_bl.txt
    black_list="pesudo_bl.txt"
    genome_size=2870182909
	promoter_file="/home/ryan/benchmark/singularity/Genome/rn6/promoter_region.bed"
	coding_promoter="/home/ryan/benchmark/singularity/Genome/rn6/coding_promoter_region.bed"
	macs2_genome="mm"
elif [[ $species == personalize ]];
	then
	echo "please add all your preferred file as reference, please make sure that you are very clear of which file is for which"
	echo "remove exit 1 after adding your file"
	exit 1
	baw_ref=" "
	chrom_size=" "
	black_list=" "
	genome_size=
	promoter_file=" "
	macs2_genome=" "
	coding_promoter=" "
fi
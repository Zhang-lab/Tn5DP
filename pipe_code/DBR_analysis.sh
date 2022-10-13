#!/bin/bash
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# group should be the character that distinguish 2 bed/peak file
group1_key=$1
group2_key=$2
# we support all software compacted in CA pipeline
cpm_value=$3
qvalue=$4
log2FC=$5

if [ -z "$group1_key" ] 
    then
    echo "missing group1 input!!!"
    echo "group_key should be the character that distinguish 2 bed/peak file"
    exit
elif [ -z "$group2_key" ] 
    then
    echo "missing group2 input!!!"
    echo "group_key should be the character that distinguish 2 bed/peak file"
    exit
elif [ -z "$cpm_value" ] || [ -z "$qvalue" ] || [ -z "$log2FC" ]
    then
    echo "missing differential analysis parameters!!!"
    echo "please select parameter values for CPM, q-value(FDR), and Log2FC!!!"
    exit
fi

# make a input count table
cat step4.1* | cut -f 1-3 > tmp_all_peak.bed
sort -k1,1V -k2,2n tmp_all_peak.bed > sorted_all_peak.bed
bedtools merge -i sorted_all_peak.bed > bedmerge_all_peak.bed

rm tmp_all_peak.bed sorted_all_peak.bed
cp bedmerge_all_peak.bed merge_counts

for file in $(ls *bed | grep $group1_key )
do
  echo $file
  intersectBed -a merge_counts -b $file -F 0.5 -wa -c > temp1.bed
  mv temp1.bed merge_counts
  echo $file >> name.txt
done

for file in $(ls *bed | grep $group2_key )
do
  echo $file
  intersectBed -a merge_counts -b $file -F 0.5 -wa -c > temp2.bed
  mv temp2.bed merge_counts
  echo $file >>name.txt
done

### step1: create the countData table 
cat bedmerge_all_peak.bed | awk -v OFS="," '{print $1,$2,$3}' > first_col_name.txt
cut -f 4- merge_counts > temp.txt
awk '{print($3-$2)}' bedmerge_all_peak.bed > temp_peak_length.txt

paste first_col_name.txt temp_peak_length.txt > temp2.txt
# create DESeq2 input count table
paste temp2.txt temp.txt > DESeq2_input_table.txt
# create edgeR input count table
cut -f 1,3- DESeq2_input_table.txt > edgeR_input_table.txt

# remove unneeded files
rm first_col_name.txt temp.txt temp_peak_length.txt temp2.txt DESeq2_input_table.txt bedmerge_all_peak.bed merge_counts name.txt

# count the sample size of each conditions
num1=$(ls *bed | grep $group1_key | wc -l)
num2=$(ls *bed | grep $group2_key | wc -l)

### step2: normalize the table by TMM method
### step3: Differential Binding Sites analysis
Rscript ${pipe_path}/DBR_analysis.R $group1_key $group2_key $num1 $num2 $cpm_value $qvalue $log2FC

















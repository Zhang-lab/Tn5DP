#!/bin/bash

# pipe start
###################################################################################################
# read all necessary parameters and prepare data structure
date
pipe_version="test_v1"
host="zhanglab/benchmark CUT&Tag-seq ATAC-seq"

# get the absolute path
pipe_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# verify data integrity using MD5 (Message Digest Algorithm 5)
md5=$(md5sum $0 | awk '{print $1}')

##################################################################################################

help_page() {
    echo "
Analyzing ATAC-seq or CUT&Tag-seq data

Program: Tn5-DP (CUT&Tag-seq and ATAC-seq analysis pipeline)
Version: 1.1
Contact: Siyuan Cheng <siyuancheng@wustl.edu>, Bo Zhang <bzhang29@wustl.edu>

Usage: path-to-pipe/pipe.sh -d  <ATAC/CUTnTag> -g <hg38/hg19/mm9/mm10/mm39/danRer11/danRer10/dm6/rn6> 
                            -r  <PE/SE>                     -m <narrow/broad>
                            -o  <experimental read file1>   -O <experimental read file2> 
                            -i  <IgG_control read file1>    -I <IgG_control read file2> 
                            -s  <AIAP/MACS2/F-seq/SEACR/HMMRATAC/ChromHMM> 
                            --personalize   (use this option when using customized peakcalling parameter setting, default is False)
                            --peakcalling options (if you use --personalize option)

Default options:        -d      data types. ATAC, CUTnTag.
                        -g      input species. Please notice that each docker/singularity image is designed for one species only.
                        -r      reads type. PE for paired-end reads, SE for single-end reads.
                        -m      marker type: narrow, broad. Please notice ATAC-seq data only support narrow marker.
                        -o      experimental read file1.
                        -O      experimental read file2.
                        -c      methylQA read length cutoff (default: 38)
                        -t      threads used for the pipe, mainly involved in cutadapt and bwa mem steps (default: 12)
                        -i      IgG control read file1.
                        -I      IgG control read file2.
                        -s      software for peakcalling: AIAP, MACS2, F-seq, SEACR, HMMRATAC, ChromHMM. default: ATAC--AIAP, CUTnTag_narrow--MACS2, CUTnTag_broad--ChromHMM
                        (Notice:        for ATAC-seq data, you could only select from AIAP(default), F-seq and HMMRATAC(only support Pair-end data) as peakcalling software,
                                        for CUT&Tag-seq data, we recommend(default) using MACS2 for narrow peakcalling, ChromHMM for broad peakcalling)
                        --personalize   users-defined parameters during peak calling step (default: False)
                        (Notice:        if user select personalized option, user need to define ALL the optional parameters for peakcalling software, 
                                        or user need set parameters as SOFTWARE-DEFINED default setting)

Personalized Peakcalling options: 

MACS2 options:          --keepdup       --keep-dup <integer> (CA default: ATAC--1000, CUT&Tag--1; MACS2-default: 1)
                        --nomodel       --nomodel (CA default: ATAC--True, CUT&Tag--False; MACS2-default: False)
                        --shift         --shift <integer> (CA default: ATAC--0, CUT&Tag--0; MACS2-default: 0)
                        --extsize       --extsize <integer> (CA default: ATAC--150, CUT&Tag--200; MACS2-default: 200)
                        --qvalue        -q | --qvalue <0~1> (CA default: ATAC--0.01, CUT&Tag--0.01; MACS2-default: 0.05)
                        --callsummits   --call-summits (CA default: False; MACS2-default: False)
                        --broad         --broad (CA default: Narrow--False, Broad--True; MACS2-default: False)
                        (Notice:        --call-summits and --broad could NOT be used at the same time)

ChromHMM options:       --binsize               -b binsize <integer> (CA default: ATAC--200, CUT&Tag--200; ChromHMM-default: 200)
                        --seed                  -s seed <integer> (CA default: ATAC--999, CUT&Tag--999; ChromHMM-default: 999)
                        --foldthresh            -f foldthresh <0~1> (CA default: ATAC--0, CUT&Tag--0; ChromHMM-default: 0)
                        --zerotransitionpower   -z zerotransitionpower <integer> (CA default: ATAC--8, CUT&Tag--8; ChromHMM-default: 8)
                        --init                  -init <information|random|load> (CA default: ATAC--information, CUT&Tag--information; ChromHMM-default: information)

HMMRATAC options        --upper         -u <integer> (CA default: ATAC--35, CUT&Tag--35; HMMRATAC-default: 20)
                        --lower         -l <integer> (CA default: ATAC--2, CUT&Tag--2; HMMRATAC-default: 10)
                        --blacklist     -e <blacklist BED file> (CA default: ATAC--none, CUT&Tag--none; HMMRATAC-default: none)
                        --score         --score <max|ave|med|fc|zscore|all> (CA default: ATAC--max, CUT&Tag--max; HMMRATAC-default: max)
                        --bedgraph      --bedgraph <True|False> (CA default: ATAC--True, CUT&Tag--True; HMMRATAC-default: False)

SEACR options           --normalize     normalize step <norm|non> (CA default: CUT&Tag--norm)
                        --model         peakcalling step <relaxed|stringent> (CA default: CUT&Tag--stringent)
"
}

# read parameters
ARGS=$(getopt -a -o d:t:g:m:o:O:r:c:i:I:s:h \
--long \
personalize,\
keepdup:,nomodel,shift:,extsize:,qvalue:,broad,callsummits,\
binsize:,seed:,foldthresh:,zerotransitionpower:,init:,\
upper:,lower:,blacklist:,score:,bedgraph,\
normalize:,model:,help \
-- "$@")

if [ $? != 0 ]; then
    echo "Terminating..." >&2
    exit 1
fi

# rearrange the order of parameters
eval set -- "${ARGS}"
# process parameters in loop by shift and while commands

main() {

    # default CA pipeline options
    threads=12
    methylQA_cutoff=38
    # swith for personalized peakcalling options
    personalize=False
    # personalized MACS2 options (default setting)
    keepdup=1
    nomodel=False
    shift=0
    extsize=200
    qvalue=0.05
    callsummits=False
    broad=False
    # personalized ChromHMM options (default setting)
    binsize=200
    seed=999
    foldthresh=0
    zerotransitionpower=8
    init=information
    # personalized HMMRATAC options (default setting)
    upper=20
    lower=10
    blacklist=
    score=max
    bedgraph=False
    # personalized SEACR options (default setting)
    normalize=
    model=

    while true; do
        case "$1" in

        # required input of default CA pipeline
        -d)
            data="$2"
            echo "assay type is: $data"
            shift 2
            ;;
        -t)
            threads="$2"
            echo "threads number for bwa is: $threads"
            shift 2
            ;;
        -g)
            species="$2"
            echo "specified species is: $species"
            shift 2
            ;;
        -o)
            R1="$2"
            echo "experimental read1: $R1"
            shift 2
            ;;
        -O)
            R2="$2"
            echo "experimental read1: $R2"
            shift 2
            ;;
        -r)
            types="$2"
            echo "type of reads is: $types"
            shift 2
            ;;
        -c)
            methylQA_cutoff="$2"
            echo "methylQA cutoff value is: $methylQA_cutoff"
            shift 2
            ;;
        -i)
            C1="$2"
            echo "IgG control read1: $C1"
            shift 2
            ;;
        -I)
            C2="$2"
            echo "IgG control read2: $C2"
            shift 2
            ;;
        -s)
            software="$2"
            echo "peakcalling software: $software"
            shift 2
            ;;
        -m)
            marker="$2"
            echo "marker type is: $marker"
            shift 2
            ;;
        # decide whether use personalized parameters or not
        --personalize)
            personalize=True
            echo "personalized peakcalling model is $personalize"
            shift
            ;;
        # MACS2 peakcalling optional choices
        --nomodel)
            nomodel=True
            echo "MACS2 --nomodel True"
            shift
            ;;
        --keepdup)
            keepdup="$2"
            echo "MACS2 --keep-dup $keepdup"
            shift 2
            ;;
        --shift)
            shift="$2"
            echo "MACS2 --shift $shift"
            shift 2
            ;;
        --extsize)
            extsize="$2"
            echo "MACS2 --extsize $extsize"
            shift 2
            ;;
        --qvalue)
            qvalue="$2"
            echo "MACS2 --qvalue $qvalue"
            shift 2
            ;;
        --callsummits)
            callsummits=True
            echo "MACS2 --call-summits True"
            shift
            ;;
        --broad)
            broad=True
            echo "MACS2 --broad True"
            shift
            ;;
        # ChromHMM peakcalling optional choices
        --foldthresh)
            foldthresh="$2"
            echo "ChromHMM -f $foldthresh"
            shift 2
            ;;
        --binsize)
            binsize="$2"
            echo "ChromHMM -b $binsize"
            shift 2
            ;;
        --seed)
            seed="$2"
            echo "ChromHMM -s $seed"
            shift 2
            ;;
        --zerotransitionpower)
            zerotransitionpower="$2"
            echo "ChromHMM -z $zerotransitionpower"
            shift 2
            ;;
        --init)
            init="$2"
            echo "ChromHMM -init $init"
            shift 2
            ;;
        # HMMRATAC peakcalling optional choices
        --upper)
            upper="$2"
            echo "HMMRATAC -u $upper"
            shift 2
            ;;
        --lower)
            lower="$2"
            echo "HMMRATAC -l $lower"
            shift 2
            ;;
        --blacklist)
            blacklist="$2"
            echo "HMMRATAC -e $blacklist"
            shift 2
            ;;
        --score)
            score="$2"
            echo "HMMRATAC --score $score"
            shift 2
            ;;
        --bedgraph)
            bedgraph=True
            echo "HMMRATAC --bedgraph $bedgraph"
            shift
            ;;
        # SEACR peakcalling optional chioces
        --normalize)
            normalize="$2"
            echo "SEACR --normalize $normalize"
            shift 2
            ;;
        --model)
            model="$2"
            echo "SEACR --model $model"
            shift 2
            ;;
        # help page and other staff...
        -h | --help)
            help_page
            exit
            ;;
        --)
            shift
            break
            ;;
        *) break ;;
        esac
    done
}

main "$@"


############################################################################################
if [ -z "$threads" ]; then
    threads=12
fi

if [ -z "$methylQA_cutoff" ]; then
    methylQA_cutoff=38
fi

if [ -z "$software" ] && [[ $marker == narrow ]] && [[ $data == CUTnTag ]]; then
    software=MACS2
elif [ -z "$software" ] && [[ $marker == broad ]] && [[ $data == CUTnTag ]]; then
    software=ChromHMM
elif [ -z "$software" ] && [[ $marker == narrow ]] && [[ $data == ATAC ]]; then
    software=AIAP
elif [[ $marker == broad ]] && [[ $data == ATAC ]]; then
    echo 'ATAC-seq does not have broad model, please use narrow model for ATAC-seq data analysis, thanks!'
    exit
fi

if [[ $R1 == *.fastq* || $R1 == *.fq.gz ]] && [[ $types == PE ]]; then
    name=$(echo ${R1%.fastq*})
    name=$(echo ${name%.fq.gz})
    raw1=$R1
    raw2=$R2
elif [[ $R1 == *.fastq* || $R1 == *.fq.gz ]] && [[ $types == SE ]]; then
    name=$(echo ${R1%.fastq*})
    name=$(echo ${name%.fq.gz})
    raw1=$R1
else
    echo "please use fastq(or fastq.gz) file......"
    exit
fi

if [[ $C1 == *.fastq* || $C1 == *.fq.gz ]] && [[ $types == PE ]] && [[ $data == CUTnTag ]]; then
    control_name=$(echo ${C1%.fastq*})
    control_name=$(echo ${control_name%.fq.gz})
    control_raw1=$C1
    control_raw2=$C2
elif [[ $C1 == *.fastq* || $C1 == *.fq.gz ]] && [[ $types == SE ]] && [[ $data == CUTnTag ]]; then
    control_name=$(echo ${C1%.fastq*})
    control_name=$(echo ${control_name%.fq.gz})
    control_raw1=$C1
fi

if [[ $C1 == *.fastq* || $C1 == *.fq.gz ]] && [[ $data == ATAC ]]; then
    echo "please do not use IgG for ATAC-seq data......"
    exit
fi

if [[ $data == ATAC ]]; then
    if [[ $software == SEACR ]]; then
        echo 'SEACR is not designed for ATAC-seq data analysis, please select approritate software...'
        exit
    elif [[ $software == HMMRATAC ]] && [[ $types == SE ]]; then
        echo 'HMMRATAC is designed for Pair-end ATAC-seq data, it could not run for Single-end data, please select approritate software...'
        exit
    elif [[ $software == ChromHMM ]]; then
        echo 'ChromHMM is not suitable for ATAC-seq data analysis, please select approritate software...'
        exit
    elif [[ $software == MACS2 ]]; then
        echo 'peakcalling software is changed to AIAP...'
        software=AIAP
    elif [[ $software == AIAP ]] && [[ $personalize == True ]]; then
        software=AIAP
    fi
elif [[ $data == CUTnTag ]]; then
    if [[ $software == AIAP ]]; then
        echo 'AIAP is not suitable for CUT&Tag-seq data analysis, please select approritate software...'
        exit
    elif [[ $software == HMMRATAC ]] && [[ $types == SE ]]; then
        echo 'HMMRATAC is designed for Pair-end data, it could not run for Single-end data, please select approritate software...'
        exit
    elif [[ $software == HMMRATAC ]] && [[ $types == PE ]]; then
        echo 'Note: HMMRATAC is designed for ATAC-seq Pair-end data, and it would not use IgG control for peakcalling...'
    elif [[ $software == ChromHMM ]] && [[ $marker == narrow ]]; then
        echo 'NOTE: ChromHMM is recommended for CUT&Tag-seq broad peakcalling, but analysis steps would continue...'
    elif [[ $software == MACS2 ]] && [[ $marker == broad ]]; then
        echo 'NOTE: MACS2 is recommended for CUT&Tag-seq narrow peakcalling, but analysis steps would continue...'
    fi
fi

if [[ $software == MACS2 ]] || [[ $software == AIAP ]]
    then
    if [[ $callsummits == True ]] && [[ $broad == True ]]
        then
        echo "--call-summits and --broad could NOT be used at the same time, please select proper peakcalling options"
        exit
    fi
fi


if [[ -z $species ]]; then
    echo 'Please select a reference genome!'
    exit
fi

if [[ -z $types ]]; then
    echo 'Please select your data type (Pair-end or Single-end)'
    exit
fi

# analysis code
# each '0' step means those are prerequest for the following one
# each step would assume previous steps have been processed
###################################################################################################

# step0, preparation (for both ATAC-seq and CUT&Tag-seq)
s0_prepare() {
    echo 'preparation step for later steps'

    mkdir 'Processed_'$name
    ln -s $(realpath $R1) ./'Processed_'$name/$R1
    ln -s $(realpath $R2 2>/dev/null) ./'Processed_'$name/$raw2 2>/dev/null
    ln -s $(realpath $C1 2>/dev/null) ./'Processed_'$name/$control_raw1 2>/dev/null
    ln -s $(realpath $C2 2>/dev/null) ./'Processed_'$name/$control_raw2 2>/dev/null
    cd ./'Processed_'$name/
    ### need to write source.sh file ###
    # source $pipe_path'/qc_source.sh' $species
    source $pipe_path'/qc_source.sh' "$species"
    ####################################
    mkdir 'QC_data_collection_'"$name"
    touch QC_pipe_processing.log

    # start record
    date >>QC_pipe_processing.log
    echo "Target file is $R1 $R2" >>QC_pipe_processing.log
    echo "assay type is $data" >>QC_pipe_processing.log

    if [[ $data == ATAC ]]; then
        echo 'there is no IgG control file for ATAC-seq data' >>QC_pipe_processing.log
    elif [[ $data == CUTnTag ]]; then
        echo "IgG control file is $C1 $C2" >>QC_pipe_processing.log
    fi
    echo "Specified species is $species" >>QC_pipe_processing.log
    echo "types of reads is $types" >>QC_pipe_processing.log
    echo "peak calling software is $software" >>QC_pipe_processing.log
    echo "marker type is $marker" >>QC_pipe_processing.log
    echo " " >>QC_pipe_processing.log
    if [[ $personalize == False ]]; then
        echo "personalized peakcalling model is $personalize, will run CA default pipeline" >>QC_pipe_processing.log
    elif [[ $personalize == True ]]; then
        echo "personalized peakcalling model is $personalize, please pay attention to the parameters" >>QC_pipe_processing.log
        if [[ $software == AIAP ]] || [[ $software == MACS2 ]]; then
            echo "MACS2 --nomodel $nomodel" >>QC_pipe_processing.log
            echo "MACS2 --keep-dup $keepdup" >>QC_pipe_processing.log
            echo "MACS2 --shift $shift" >>QC_pipe_processing.log
            echo "MACS2 --extsize $extsize" >>QC_pipe_processing.log
            echo "MACS2 --qvalue $qvalue" >>QC_pipe_processing.log
            echo "MACS2 --call-summits $callsummits" >>QC_pipe_processing.log
            echo "MACS2 --broad $broad" >>QC_pipe_processing.log
        elif [[ $software == ChromHMM ]]; then
            echo "ChromHMM --foldthresh $foldthresh" >>QC_pipe_processing.log
            echo "ChromHMM --binsize $binsize" >>QC_pipe_processing.log
            echo "ChromHMM --seed $seed" >>QC_pipe_processing.log
            echo "ChromHMM --zerotransitionpower $zerotransitionpower" >>QC_pipe_processing.log
            echo "ChromHMM --init $init" >>QC_pipe_processing.log
        elif [[ $software == HMMRATAC ]]; then
            echo "HMMRATAC --upper $upper" >>QC_pipe_processing.log
            echo "HMMRATAC --lower $lower" >>QC_pipe_processing.log
            echo "HMMRATAC --blacklist $blacklist" >>QC_pipe_processing.log
            echo "HMMRATAC --score $score" >>QC_pipe_processing.log
            echo "HMMRATAC --bedgraph $bedgraph" >>QC_pipe_processing.log
        elif [[ $software == SEACR ]]; then
            echo "SEACR --normalize $normalize" >>QC_pipe_processing.log
            echo "SEACR --model $model" >>QC_pipe_processing.log
        fi
    fi
    echo " " >>QC_pipe_processing.log

    if [[ $R1 == *.sra ]]; then
        $fastq_dump_tool $R1 --split-3
        if [[ $types == PE ]]; then
            raw1=$name'_1.fastq'
            raw2=$name'_2.fastq'
        elif [[ $types == SE ]]; then
            raw1=$name'.fastq'
        fi
    fi
}

# Step 1.1, Trim ATAC-seq/CUT&Tag-seq adapters and QC on seq file
s1.1_cutadapt() {
    echo 'cutadapt step1.1'

    if [[ $data == ATAC ]] && [[ $types == PE ]]; then
        echo 'trimming ATAC PE by cutadapt'
        $cutadapt -a $adapter_1 -A $adapter_2 --quality-cutoff=15,10 --minimum-length=25 -o 'step1.1_trimed_'$name'_1.fastq' -p 'step1.1_trimed_'$name'_2.fastq' $raw1 $raw2 >'step1.1_'$name'_cutadapt_PE.trimlog'
        temp=$(grep "Total read pairs processed:" step1.1_*trimlog | awk '{print $5}')
        raw_reads=$(echo ${temp//,/})
        temp2=$(grep "Pairs written" step1.1_*trimlog | awk '{print $5}')
        written_reads=$(echo ${temp2//,/})
    elif [[ $data == ATAC ]] && [[ $types == SE ]]; then
        echo 'trimming ATAC SE reads by cutadapt'
        $cutadapt -a $adapter_1 --quality-cutoff=15,10 --minimum-length=25 -o 'step1.1_trimed_'$name'.fastq' $raw1 >'step1.1_'$name'_cutadapt_SE.trimlog'
        temp=$(grep "Total reads processed:" step1.1_*trimlog | awk '{print $4}')
        raw_reads=$(echo ${temp//,/})
        temp2=$(grep "Reads written" step1.1_*trimlog | awk '{print $5}')
        written_reads=$(echo ${temp2//,/})

    elif [[ $data == CUTnTag ]] && [[ $types == PE ]]; then
        # process input raw reads
        echo 'trimming CUT&Tag input PE reads by cutadapt'
        $cutadapt -a "$adapter_1" -A "$adapter_2" --quality-cutoff=15,10 --minimum-length=25 -o 'step1.1_trimed_'$name'_1.fastq' -p 'step1.1_trimed_'$name'_2.fastq' $raw1 $raw2 >'step1.1_'$name'_cutadapt_PE.trimlog'
        temp=$(grep "Total read pairs processed:" step1.1_"$name"*trimlog | awk '{print $5}')
        raw_reads=$(echo ${temp//,/})
        temp2=$(grep "Pairs written" step1.1_"$name"*trimlog | awk '{print $5}')
        written_reads=$(echo ${temp2//,/})
        # process IgG control reads
        echo 'trimming CUT&Tag IgG_control PE reads by cutadapt'
        $cutadapt -a "$adapter_1" -A "$adapter_2" --quality-cutoff=15,10 --minimum-length=25 -o 'step1.1_trimed_'$control_name'_1.fastq' -p 'step1.1_trimed_'$control_name'_2.fastq' $control_raw1 $control_raw2 >'step1.1_'$control_name'_cutadapt_PE.trimlog'

    elif [[ $data == CUTnTag ]] && [[ $types == SE ]]; then
        # process input raw reads
        echo 'trimming CUT&Tag input SE reads by cutadapt'
        $cutadapt -a ${adapter_1} --quality-cutoff=15,10 --minimum-length=25 -o 'step1.1_trimed_'$name'.fastq' $raw1 >'step1.1_'$name'_cutadapt_SE.trimlog'
        temp=$(grep "Total reads processed:" step1.1_"$name"*trimlog | awk '{print $4}')
        raw_reads=$(echo ${temp//,/})
        temp2=$(grep "Reads written" step1.1_"$name"*trimlog | awk '{print $5}')
        written_reads=$(echo ${temp2//,/})
        # process IgG control raw reads
        echo 'trimming CUT&Tag IgG_control SE reads by cutadapt'
        ${cutadapt} -a ${adapter_1} --quality-cutoff=15,10 --minimum-length=25 -o 'step1.1_trimed_'$control_name'.fastq' $control_raw1 >'step1.1_'$control_name'_cutadapt_SE.trimlog'
    fi

    if [ $? == 0 ]; then
        echo "step1.1, cutadapt trimming done" >>QC_pipe_processing.log
    else
        echo "step1.1, cutadapt trimming fail......" >>QC_pipe_processing.log
        # exit 1
    fi
}

s1.2_fastqc() {
    echo 'fastqc is processing fastq file......'
    if [[ $data == ATAC ]]; then
        [ -f $(ls 'step1.1_trimed_'$name*'.fastq' | head -1) ] && $fastqc -t $threads 'step1.1_trimed_'$name*'.fastq' -o .
    elif [[ $data == CUTnTag ]]; then
        [ -f $(ls 'step1.1_trimed_'$name*'.fastq' | head -1) ] && $fastqc -t $threads 'step1.1_trimed_'$name*'.fastq' -o .
        [ -f $(ls 'step1.1_trimed_'$control_name*'.fastq' | head -1) ] && $fastqc -t $threads 'step1.1_trimed_'$control_name*'.fastq' -o .
    fi

    if [ $? == 0 ]; then
        echo "step1.2, fastqc process done" >>QC_pipe_processing.log
    else
        echo "step1.2, fastqc process fail......" >>QC_pipe_processing.log
        exit 1
    fi

    for zip in $(ls | grep fastqc.zip); do
        unzip -o $zip
        mv $zip 'step1.2_'$zip
    done

    # 1.3 fastqc data collection #
    # only collect the fastqc data of experimental data
    echo -e "filename\tdeduplication_percentage\tmarker" >'step1.3_dedup_percentage_'$name'.result'
    for file in $(ls -d *$name*'fastqc'/); do
        cd $file
        temp=$(echo ${file##step1.1_trimed_})
        out_name=$(echo ${temp%*_fastqc/})
        out_value=$(grep 'Total Deduplicated Percentage' fastqc_data.txt | awk '{print $4}')
        echo -e "$out_name\t$out_value\t$marker" >>../'step1.3_dedup_percentage_'$name'.result'
        echo -e "item\t$out_name\t$out_name" >'step1.3_duplication_summary_'$out_name'.result'
        grep 'Sequence Duplication Levels' -A 15 fastqc_data.txt >>'step1.3_duplication_summary_'$out_name'.result'
        mv 'step1.3_duplication_summary_'$out_name'.result' ../'QC_data_collection_'$name
        echo -e "$out_name\tfastqc_test" >'step1.3_fastqc_summary_'$out_name'.result'
        awk -F "\t" '{print $1,$2}' OFS='\t' summary.txt >>'step1.3_fastqc_summary_'$out_name'.result'
        mv 'step1.3_fastqc_summary_'$out_name'.result' ../'QC_data_collection_'$name
        cd ..
    done
    # step 1.3 data cleaning before next step
    if [ $? == 0 ]; then
        echo "step1.3, fastqc data_collection process done" >>QC_pipe_processing.log
    else
        echo "step1.3, fastqc data_collection process fail......" >>QC_pipe_processing.log
    fi

    sed 1d step1.3_dedup_percentage_$name'.result' | cut -f 2 >temp_dedup.txt &&
        before_dedup=$(python3 -c "print($(awk '{s+=$1}END{print s}' temp_dedup.txt) * 0.01 /$(cat temp_dedup.txt | wc -l))") &&
        before_dup=$(python3 -c "print(1-$before_dedup*1.0)") &&
        rm temp_dedup.txt
    mv 'step1.3_dedup_percentage_'$name'.result' ./'QC_data_collection_'$name
    mv *fastqc* 'QC_data_collection_'$name

    # 1.4, get PE data R1 R2 deduplication difference percentage
    # only collect the fastqc data of experimental data
    if [[ $types == PE ]]; then
        per1=$(tail -n 2 ./'QC_data_collection_'$name/'step1.3_dedup_percentage_'$name'.result' | awk '{print $2}' | sed -n '1p')
        per2=$(tail -n 2 ./'QC_data_collection_'$name/'step1.3_dedup_percentage_'$name'.result' | awk '{print $2}' | sed -n '2p')
        dif=$(echo "scale=2; ($per1-$per2)*200/($per1+$per2)" | bc -l)
    else
        dif=0
    fi

    if [ $? == 0 ]; then
        echo "step1.4, calculate replicate difference process done" >>QC_pipe_processing.log
    else
        echo "step1.4, calculate replicate difference process fail......" >>QC_pipe_processing.log
    fi
}

# step2.0, files check
s2.0_ref() {
    echo 'step2.0 files check'

    # refine chrom_size file (remove random and Unknown record)
    awk '{ if ((length($1) < 6) && (length($1) > 1))  print $0}' OFS='\t' $chrom_size >refined_chrom_size.txt
    chrom_size=$(pwd)"/refined_chrom_size.txt"
}

# step2.1, BWA MEM alignment
s2.1_bwa() {
    echo 'alignment by bwa......'
    if [[ $data == ATAC ]]; then
        $bwa mem -t $threads $bwa_ref 'step1.1_trimed_'$name*'.fastq' | $samtools view -bS - | $samtools sort - -O 'bam' -o 'step2.1_trimed_'$name'.bam' -T temp_aln_input
    elif [[ $data == CUTnTag ]]; then
        $bwa mem -t $threads $bwa_ref 'step1.1_trimed_'$name*'.fastq' | $samtools view -bS - | $samtools sort - -O 'bam' -o 'step2.1_trimed_'"$name"'.bam' -T temp_aln_input
        $bwa mem -t $threads $bwa_ref 'step1.1_trimed_'$control_name*'.fastq' | $samtools view -bS - | $samtools sort - -O 'bam' -o 'step2.1_trimed_'$control_name'.bam' -T temp_aln_IgG
    fi

    if [ $? == 0 ]; then
        echo "step2.1, bwa alignment process done" >>QC_pipe_processing.log
        rm 'step1.1_trimed_'$name*'.fastq'
    else
        echo "step2.1, bwa alignment process fail......" >>QC_pipe_processing.log
        exit 1
    fi
}

# step2.2, removing low mapQ reads and count reads distribution (mapQ=0 makes no sense here, because they are not reliable)
# only focus on Marker_input data, don not calculate low mapQ reads of IgG control
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

# step2.3, preseq
s2.3_preseq() {
    echo 'step2.3 preseq...'

    if [[ $data == ATAC ]]; then
        $preseq lc_extrap -o 'step2.3_yield_'$name'.result' -B 'step2.1_trimed_'$name'.bam'
        if [ $? == 0 ]; then
            echo "step2.3, preseq lc_extrap estimate process done" >>QC_pipe_processing.log
        else
            echo "step2.3, preseq lc_extrap estimate process fail......" >>QC_pipe_processing.log
        fi
        mv 'step2.3_yield_'$name'.result' ./'QC_data_collection_'$name

    elif [[ $data == CUTnTag ]]; then
        $preseq lc_extrap -o 'step2.3_yield_'$name'.result' -B 'step2.1_trimed_'$name'.bam'
        $preseq lc_extrap -o 'step2.3_yield_'$control_name'.result' -B 'step2.1_trimed_'$control_name'.bam'
        if [ $? == 0 ]; then
            echo "step2.3, preseq lc_extrap estimate process done" >>QC_pipe_processing.log
        else
            echo "step2.3, preseq lc_extrap estimate process fail......" >>QC_pipe_processing.log
        fi
        mv 'step2.3_yield_'$name'.result' ./'QC_data_collection_'$name
        mv 'step2.3_yield_'$control_name'.result' ./'QC_data_collection_'$name
    fi
}

# 3.1, methylQA
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

# run pipe ()
###################################################################################################
# step-by-step
s0_prepare
s1.1_cutadapt
s1.2_fastqc
s2.0_ref
s2.1_bwa
s2.2_distri
s2.3_preseq
s3.1_methylQA
s3.2_rmbl_bg
s4.0_IgG_data_pre
s4.1_peakcall

echo "Processing $name done"
cd ..
date

FROM ubuntu:16.04

### AC base for v1.0

### excutables that CAN NOT be controlled
# 2. methylQA: no version control for install file (07/06, v0.2.1)

### freeze tool version
# Python 3.9.6
# cutadapt 3.5
# R 4.1.2
# macs2 2.2.7
# Bioconductor 3.14 (BiocManager 1.30.18)
# SRA toolkit (fastq-dump) 2.10.0
# fastqc 0.11.7
# preseq 2.0.0
# bwa 0.7.16a
# samtools 1.9
# bedtools 2.29
# methylQA: no version freeze
## software not installed in AIAP
# Java 1.8
# perl 5.22.1
# F-seq 1.85
# SEACR 1.3
# HMMRATAC 1.2.10
# ChromHMM 1.23


# install all build-essentials
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        bc \
        bzip2 \
        ca-certificates \
        curl \
        g++ \
        gcc \
        gfortran \
        git \
        less \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgsl-dev \
        libgsl2 \
        liblzma-dev \
        libmysqlclient-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        libx11-dev \
        libxml2 \
        libxml2-dev \
        make \
        parallel \
        python-dev \
        texlive-latex-base texlive-latex-extra \
        unzip \
        xorg-dev \
        zlib1g-dev \
        vim \
        wget \
        python-pip \
        python-setuptools \
    && rm -rf /var/lib/apt/lists/* 


# install python3.9.6 and cutadapt 3.5 
RUN wget https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz --no-check-certificate -o /opt/Python-3.6.3.tgz \
    && cd /opt && tar -xvf Python-3.9.6.tgz \
    && cd Python-3.9.6 \
    && ./configure \
    && make \
    && make install \
    && rm /opt/Python-3.9.6.tgz \
    && apt-get update \
    && apt-get install -y python3-pip \
    && pip install --upgrade pip \
    && pip install 'cutadapt==3.5' \
    && pip install  numpy scipy multiqc \
    && pip install 'MACS2==2.2.7'


# install R 4.1.2
RUN cd /opt/ \
    && wget -c https://repo.continuum.io/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh --no-check-certificate -o /opt/Miniconda3-py39_4.12.0-Linux-x86_64.sh\
    && bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -p /opt/miniconda3 \
    && rm Miniconda3-py39_4.12.0-Linux-x86_64.sh \
    && ln -s /opt/miniconda3/bin/conda /usr/bin/ \
    && conda config --add channels defaults \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda install -y -c r r-base=4.1.2 

ENV PATH $PATH:/opt/miniconda3/bin/
### to test: 
# export PATH=$PATH:/opt/miniconda3/bin/

RUN echo 'install.packages("ggplot2",repos="http://cran.us.r-project.org")' > /opt/packages.R \
    && echo 'install.packages("cowplot",repos="http://cran.us.r-project.org")' >> /opt/packages.R \
    && echo 'install.packages("jsonlite",repos="http://cran.us.r-project.org")' >> /opt/packages.R \
    && echo 'install.packages("data.table",repos="http://cran.us.r-project.org")' >> /opt/packages.R \
    && echo 'install.packages("BiocManager",repos="http://cran.us.r-project.org")' >> /opt/packages.R \
    && echo 'BiocManager::install(version = "3.14")' >> /opt/packages.R \
    && echo 'BiocManager::install("DESeq2")' >> /opt/packages.R \
    && echo 'BiocManager::install("edgeR")' >> /opt/packages.R \
    && echo 'BiocManager::install("DiffBind")' >> /opt/packages.R
    
RUN /opt/miniconda3/bin/Rscript /opt/packages.R \
    && rm /opt/packages.R \
    && mkdir -p /usr/local/lib/R/site-library

# install sra toolkit 2.10.0 (fastq-dump.2.10.0)
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz --no-check-certificate -o /opt/sratoolkit.current-ubuntu64.tar.gz \
    && tar -xvzf /opt/sratoolkit.current-ubuntu64.tar.gz -C /opt/ \
    && rm /opt/sratoolkit.current-ubuntu64.tar.gz

# install fastqc
RUN conda install -y -c bioconda fastqc=0.11.7

# install preseq 2.0.0
RUN wget http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 --no-check-certificate -o /opt/preseq_linux_v2.0.tar.bz2 \
    && tar xvjf /opt/preseq_linux_v2.0.tar.bz2 -C /opt/ \
    && ln -s /opt/preseq_v2.0/preseq /usr/local/bin/preseq \
    && ln -s /opt/preseq_v2.0/bam2mr /usr/local/bin/bam2mr \
    && rm /opt/preseq_linux_v2.0.tar.bz2 \
    && ln -s /usr/lib/x86_64-linux-gnu/libgsl.so /lib/x86_64-linux-gnu/libgsl.so.0
    # make sure that libgsl.so.0 exists beause PreSeq links to that

# install bwa 0.7.16a
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.16a.tar.bz2 --no-check-certificate -o /opt/bwa-0.7.16a.tar.bz2 \
    && tar xvjf /opt/bwa-0.7.16a.tar.bz2 -C /opt/ \
        && cd /opt/bwa-0.7.16a/ \
        && make \
        && rm /opt/bwa-0.7.16a.tar.bz2

# install samtools v1.9
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 -o /opt/samtools-1.9.tar.bz2 \
    && tar xvjf /opt/samtools-1.9.tar.bz2 -C /opt/ \
    && cd /opt/samtools-1.9 \
    && make \
    && make install \
    && rm /opt/samtools-1.9.tar.bz2

# install bedtools 2.29.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz -o /opt/bedtools-2.29.0.tar.gz \
        && tar -zxvf /opt/bedtools-2.29.0.tar.gz -C /opt/ \
        && cd /opt/bedtools2 \
        && make \
        && rm /opt/bedtools-2.29.0.tar.gz

# install methylQA 0.2.1
RUN cd /opt/ \
        && git clone git://github.com/lidaof/methylQA.git \
        && cd /opt/methylQA/ \
        && make \
        && cd /

# install bedGraphToBigWig
RUN conda install -y ucsc-bedgraphtobigwig

# install jq for image finding
RUN conda install -y -c conda-forge jq


# software not installed in AIAP
# install Java 1.8
RUN wget http://regmedsrv1.wustl.edu/Public_SPACE/ryan/Public_html/singularity_ac/jdk-8u131-linux-x64.tar.gz  --no-check-certificate -o /opt/jdk-18_linux-x64_bin.tar.gz \
    && cd /opt && tar -xvf jdk-18_linux-x64_bin.tar.gz \
    && rm jdk-18_linux-x64_bin.tar.gz \
    && ln -s /opt/jdk1.8.0_131/bin/java /usr/local/bin/java


# install Perl 5.22.1
RUN wget https://www.cpan.org/src/5.0/perl-5.22.1.tar.gz  --no-check-certificate -o /opt/perl-5.22.1.tar.gz \
    && cd /opt && tar -xvf perl-5.22.1.tar.gz \
    && rm perl-5.22.1.tar.gz \
    && cd  perl-5.22.1 \
    && ./Configure -de \
    && make \
    && make install


# install F-seq 1.85
RUN cd /opt/ \
    && git clone https://github.com/aboyle/F-seq.git \
    && cd F-seq \
    && ant \
    && cd dist~ \
    && tar -xvf fseq.tgz \
    && ln -s /opt/F-seq/dist~/fseq/bin/fseq /usr/local/bin/fseq


# install HMMRATAC 1.2.10
RUN cd /opt/ \
    && git clone https://github.com/LiuLabUB/HMMRATAC.git \
    && cd HMMRATAC \
    && wget https://github.com/LiuLabUB/HMMRATAC/releases/download/1.2.10/HMMRATAC_V1.2.10_exe.jar --no-check-certificates


# install SEACR 1.3
RUN cd /opt/ \
    && git clone https://github.com/FredHutch/SEACR.git


# install ChromHMM 1.23
RUN cd /opt/ \
    && wget https://ernstlab.biolchem.ucla.edu/ChromHMM/ChromHMM.zip \
    && unzip ChromHMM.zip 


ENV PATH $PATH:/opt/sratoolkit.2.10.0-ubuntu64/bin/:/opt/bwa-0.7.16a/:/opt/bedtools2/bin/:/opt/methylQA/:/opt/kentUtils/bin/













FROM ubuntu:18.04
MAINTAINER bhaas@broadinstitute.org

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
    fort77 \
    ftp \
    g++ \
    gcc \
    gfortran \
    git \
    libblas-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libghc-zlib-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    nano \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-distutils \
    python3-pip \
    python3-setuptools \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev


## Perl stuff
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install DB_File
RUN cpanm install URI::Escape


## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin


#####
# Install R

WORKDIR $SRC

ENV R_VERSION=R-3.6.3

RUN curl https://cran.r-project.org/src/base/R-3/$R_VERSION.tar.gz -o $R_VERSION.tar.gz && \
    tar xvf $R_VERSION.tar.gz && \
    cd $R_VERSION && \
	./configure && make && make install

RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("tidyverse")'
RUN R -e 'BiocManager::install("edgeR")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("ape")'
RUN R -e 'BiocManager::install("ctc")'
RUN R -e 'BiocManager::install("gplots")'
RUN R -e 'BiocManager::install("Biobase")'
RUN R -e 'BiocManager::install("qvalue")'
RUN R -e 'BiocManager::install("goseq")'
RUN R -e 'BiocManager::install("Glimma")'
RUN R -e 'BiocManager::install("ROTS")'
RUN R -e 'BiocManager::install("GOplot")'
RUN R -e 'BiocManager::install("argparse")'
RUN R -e 'BiocManager::install("fastcluster")'
RUN R -e 'BiocManager::install("DEXSeq")'
RUN R -e 'BiocManager::install("tximport")'
RUN R -e 'BiocManager::install("tximportData")'


##RUN R -e 'BiocManager::install("sm")'
##RUN R -e 'BiocManager::install("vioplot")'


ENV LD_LIBRARY_PATH=/usr/local/lib

RUN apt-get install -y cython3

## Python 3 stuff
RUN ln -sf /usr/bin/python3 /usr/bin/python

## some python modules
RUN pip3 install numpy
RUN pip3 install git+https://github.com/ewels/MultiQC.git
RUN pip3 install HTSeq


## bowtie
WORKDIR $SRC
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.1.1/bowtie-1.2.1.1-linux-x86_64.zip/download -O bowtie-1.2.1.1-linux-x86_64.zip && \
    unzip bowtie-1.2.1.1-linux-x86_64.zip && \
	mv bowtie-1.2.1.1/bowtie* $BIN

## RSEM
RUN mkdir /usr/local/lib/site_perl
WORKDIR $SRC
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz && \
    tar xvf v1.3.3.tar.gz && \
    cd RSEM-1.3.3 && \
    make && \
    cp rsem-* convert-sam-for-rsem $BIN && \
    cp rsem_perl_utils.pm /usr/local/lib/site_perl/ && \
    cd ../ && rm -r RSEM-1.3.3



## Kallisto
WORKDIR $SRC
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar xvf kallisto_linux-v0.46.1.tar.gz && \
    mv kallisto/kallisto $BIN


## FASTQC
WORKDIR $SRC
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    chmod 755 /usr/local/src/FastQC/fastqc && \
    ln -s /usr/local/src/FastQC/fastqc $BIN/.

# blast
WORKDIR $SRC
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/ncbi-blast-2.5.0+-x64-linux.tar.gz && \
    tar xvf ncbi-blast-2.5.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.5.0+/bin/* $BIN && \
    rm -r ncbi-blast-2.5.0+

## Bowtie2
WORKDIR $SRC
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip/download -O bowtie2-2.3.4.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.4.1-linux-x86_64.zip && \
    mv bowtie2-2.3.4.1-linux-x86_64/bowtie2* $BIN && \
    rm *.zip && \
    rm -r bowtie2-2.3.4.1-linux-x86_64

## Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar xvf samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && make && make install


## Jellyfish
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz && \
    tar xvf jellyfish-2.2.7.tar.gz && \
    cd jellyfish-2.2.7/ && \
    ./configure && make && make install

## Salmon
WORKDIR $SRC
ENV SALMON_VERSION=1.5.2
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    tar xvf Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    ln -s $SRC/Salmon-latest_linux_x86_64/bin/salmon $BIN/.


## FeatureCounts
RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz/download -O subread-2.0.0-Linux-x86_64.tar.gz && \
    tar xvf subread-2.0.0-Linux-x86_64.tar.gz && \
    cp -r subread-2.0.0-Linux-x86_64/bin/* $BIN/

## Hisat2
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
    unzip hisat2-2.1.0-Linux_x86_64.zip && \
    cp hisat2-2.1.0/hisat2* $BIN/

## GMAP
ENV GSNAP_VER 2021-07-23
WORKDIR $SRC
RUN GMAP_URL="http://research-pub.gene.com/gmap/src/gmap-gsnap-$GSNAP_VER.tar.gz" && \
    wget $GMAP_URL && \
    tar xvf gmap-gsnap-$GSNAP_VER.tar.gz && \
    cd gmap-$GSNAP_VER && ./configure --prefix=`pwd` && make && make install && \
    cp bin/* $BIN/


# patch for salmon
WORKDIR $SRC
RUN ln -sf /usr/local/src/salmon-latest_linux_x86_64/bin/salmon $BIN/.

# blat
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat -P $BIN && \
    chmod a+x $BIN/blat

## Picard tools
WORKDIR $SRC
RUN wget https://github.com/broadinstitute/picard/releases/download/2.23.3/picard.jar
ENV PICARD_HOME $SRC

####
## GATK4 installation
WORKDIR $SRC
ENV GATK_VERSION=4.1.8.0
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip

ENV GATK_HOME $SRC/gatk-${GATK_VERSION}


## STAR
ENV STAR_VERSION=2.7.8a
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
        tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
            mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin



## Trinity
WORKDIR $SRC
ENV TRINITY_VERSION="2.12.0"
ENV TRINITY_CO="a1431ffa1d87236dba814ddf0bb2e5c1f45c43d2"

WORKDIR $SRC

RUN git clone --recursive https://github.com/trinityrnaseq/trinityrnaseq.git && \
    cd trinityrnaseq && \
    git checkout ${TRINITY_CO} && \
    git submodule init && git submodule update && \
    git submodule foreach --recursive git submodule init && \
    git submodule foreach --recursive git submodule update && \
    rm -rf ./trinity_ext_sample_data && \
    make && make plugins && \
    make install && \
    cd ../ && rm -r trinityrnaseq

ENV TRINITY_HOME /usr/local/bin

ENV PATH=${TRINITY_HOME}:${PATH}


# some cleanup
WORKDIR $SRC
RUN rm -r ${R_VERSION} *.tar.gz *.zip *.bz2

RUN apt-get -qq -y remove git && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

################################################
## be sure this is last!
COPY Dockerfile $SRC/Dockerfile.$TRINITY_VERSION


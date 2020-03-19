FROM ubuntu:18.04
MAINTAINER bhaas@broadinstitute.org

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y gcc g++ perl python3 automake make \
                                       wget git curl libdb-dev \
                                       zlib1g-dev bzip2 libncurses5-dev \
				       texlive-latex-base \
                                       default-jre \
				       python3-pip python3-dev python3-setuptools python3-distutils \
				       gfortran \
				       build-essential libghc-zlib-dev libncurses-dev libbz2-dev liblzma-dev libpcre3-dev libxml2-dev \
				       libblas-dev gfortran git unzip ftp libzmq3-dev nano ftp fort77 libreadline-dev \
				       libcurl4-openssl-dev libx11-dev libxt-dev \
				       x11-common libcairo2-dev libpng-dev libreadline-dev libjpeg-dev pkg-config libtbb-dev \
					   cmake rsync libssl-dev tzdata \
	&& apt-get clean




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
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz && \
    tar xvf v1.3.0.tar.gz && \
    cd RSEM-1.3.0 && \
    make && \
    cp rsem-* $BIN && \
    cp rsem_perl_utils.pm /usr/local/lib/site_perl/ && \
    cd ../ && rm -r RSEM-1.3.0


## Kallisto
WORKDIR $SRC
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz && \
    tar xvf kallisto_linux-v0.43.1.tar.gz && \
    mv kallisto_linux-v0.43.1/kallisto $BIN


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


## RSEM
WORKDIR $SRC
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz && \
	 tar xvf v1.3.0.tar.gz && \
     cd RSEM-1.3.0 && \
     make && \
     cp rsem-* $BIN && \
     cp convert-sam-for-rsem $BIN && \
     cp rsem_perl_utils.pm /usr/local/lib/site_perl/ && \
     cd ../ && rm -r RSEM-1.3.0


## Salmon
WORKDIR $SRC
ENV SALMON_VERSION=1.0.0
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    tar xvf Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    ln -s $SRC/Salmon-latest_linux_x86_64/bin/salmon $BIN/.
    



##############
## STAR

ENV STAR_VERSION=2.7.2b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
        tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
            mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin

#####
## FeatureCounts

RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz/download -O subread-2.0.0-Linux-x86_64.tar.gz && \
    tar xvf subread-2.0.0-Linux-x86_64.tar.gz && \
    cp -r subread-2.0.0-Linux-x86_64/bin/* $BIN/



#####
## Hisat2

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
    unzip hisat2-2.1.0-Linux_x86_64.zip && \
    cp hisat2-2.1.0/hisat2* $BIN/


####
## GMAP

ENV GMAP_VERSION=2017-11-15
WORKDIR $SRC
RUN GMAP_URL="http://research-pub.gene.com/gmap/src/gmap-gsnap-$GMAP_VERSION.tar.gz" && \
    wget $GMAP_URL && \
    tar xvf gmap-gsnap-$GMAP_VERSION.tar.gz && \
    cd gmap-$GMAP_VERSION && ./configure --prefix=`pwd` && make && make install


####
## Picard tools

WORKDIR $SRC
RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.3/picard.jar
ENV PICARD_HOME $SRC

####
## GATK4 installation

WORKDIR $SRC

ENV GATK_VERSION=4.1.4.0
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip

ENV GATK_HOME $SRC/gatk-${GATK_VERSION}





# patch for salmon
WORKDIR $SRC
RUN ln -sf /usr/local/src/salmon-latest_linux_x86_64/bin/salmon $BIN/.

# blat
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat -P $BIN && \
    chmod a+x $BIN/blat


##########
## Trinity

WORKDIR $SRC

ENV TRINITY_VERSION="2.10.0"
ENV TRINITY_CO="3b7f0dff1ff5bbe5de6c64c7df2f78b98bcd21f1"

WORKDIR $SRC

RUN git clone --recursive https://github.com/trinityrnaseq/trinityrnaseq.git && \
    cd trinityrnaseq && \
    git checkout ${TRINITY_CO} && \
    git submodule init && git submodule update && \
    rm -rf ./trinity_ext_sample_data && \
    make && make plugins && \
    make install && \
    cd ../ && rm -r trinityrnaseq

ENV TRINITY_HOME /usr/local/bin/trinityrnaseq

ENV PATH=${TRINITY_HOME}:${PATH}





# some cleanup
WORKDIR $SRC
RUN rm -r ${R_VERSION} *.tar.gz *.zip *.bz2

RUN apt-get clean

################################################
## be sure this is last!
COPY Dockerfile $SRC/Dockerfile.$TRINITY_VERSION  

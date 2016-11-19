FROM ubuntu:16.04
MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget git curl libdb-dev \
                                       zlib1g-dev bzip2 libncurses5-dev \
				       texlive-latex-base \
                                       default-jre \
				       python-pip python-dev \
				       gfortran \
				       build-essential libghc-zlib-dev libncurses-dev libbz2-dev liblzma-dev libpcre3-dev libxml2-dev \
				       libblas-dev gfortran git unzip ftp libzmq3-dev nano ftp fort77 libreadline-dev \
				       libcurl4-openssl-dev libx11-dev libxt-dev \
				       x11-common libcairo2-dev libpng12-dev libreadline6-dev libjpeg8-dev pkg-config \
                   && apt-get clean

RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN cpanm install DB_File

RUN cpanm install URI::Escape


## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin


#####
# Install R

WORKDIR $SRC

RUN curl https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz -o R-3.2.1.tar.gz && \
        tar xvf R-3.2.1.tar.gz && \
        cd R-3.2.1 && \
	./configure && make && make install
	    

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");library(BiocInstaller); biocLite("edgeR", dep = TRUE)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");library(BiocInstaller); biocLite("DESeq2", dep = TRUE)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");library(BiocInstaller); biocLite("ape", dep = TRUE)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");library(BiocInstaller); biocLite("ctc", dep = TRUE)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");library(BiocInstaller); biocLite("gplots", dep = TRUE)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");library(BiocInstaller); biocLite("Biobase", dep = TRUE)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");library(BiocInstaller); biocLite("qvalue", dep = TRUE)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");library(BiocInstaller); biocLite("goseq", dep = TRUE)'


## bowtie
WORKDIR $SRC
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip/download -O bowtie-1.1.2-linux-x86_64.zip && \
        unzip bowtie-1.1.2-linux-x86_64.zip && \
	cp bowtie-1.1.2/bowtie* $BIN


## RSEM
RUN mkdir /usr/local/lib/site_perl
WORKDIR $SRC
RUN wget https://github.com/deweylab/RSEM/archive/v1.2.31.tar.gz && \
    tar xvf v1.2.31.tar.gz && \
    cd RSEM-1.2.31 && \
    make && \
    cp rsem-* $BIN && \
    cp rsem_perl_utils.pm /usr/local/lib/site_perl/


## Kallisto
WORKDIR $SRC
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.42.5/kallisto_linux-v0.42.5.tar.gz && \
    tar xvf kallisto_linux-v0.42.5.tar.gz && \
    cp kallisto_linux-v0.42.5/kallisto $BIN

## FASTQC
WORKDIR $SRC
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    chmod 755 /usr/local/src/FastQC/fastqc && \
    ln -s /usr/local/src/FastQC/fastqc $BIN/.


## Bowtie2
WORKDIR $SRC
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip/download -O bowtie2-2.2.9-linux-x86_64.zip && \
    unzip bowtie2-2.2.9-linux-x86_64.zip && \
    mv bowtie2-2.2.9/bowtie2* $BIN && \
    rm *.zip && \
    rm -r bowtie2-2.2.9

# blast
WORKDIR $SRC
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/ncbi-blast-2.5.0+-x64-linux.tar.gz && \
    tar xvf ncbi-blast-2.5.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.5.0+/bin/* $BIN



##########
## Trinity


WORKDIR $SRC

ENV TRINITY_CO=3953dd1
RUN git clone https://github.com/trinityrnaseq/trinityrnaseq.git && \
    cd trinityrnaseq && \
    git checkout 3953dd1
    

ENV TRINITY_HOME $SRC/trinityrnaseq

RUN cd $TRINITY_HOME && make && make plugins

RUN cp $TRINITY_HOME/trinity-plugins/htslib/bgzip $BIN

RUN cp $TRINITY_HOME/trinity-plugins/BIN/samtools $BIN

RUN cp $TRINITY_HOME/trinity-plugins/htslib/tabix $BIN

ENV PATH=${TRINITY_HOME}:${PATH}




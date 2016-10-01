FROM debian
MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget git curl libdb-dev \
                                       zlib1g-dev bzip2 libncurses5-dev \
				       texlive-latex-base \
                                       openjdk-7-jre \
				       python-pip python-dev \
    && apt-get clean

RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN cpanm install DB_File

RUN cpanm install URI::Escape

RUN pip install pysam


## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin


######################
## Tool installations:
######################

##########
## Trinity


RUN TRINITY_URL="https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz" && \
    wget -P $SRC $TRINITY_URL && \
    tar -xvf $SRC/v2.1.1.tar.gz -C $SRC && \
    cd $SRC/trinityrnaseq-2.1.1 && make

ENV TRINITY_HOME $SRC/trinityrnaseq-2.1.1


RUN cp $TRINITY_HOME/trinity-plugins/htslib/bgzip $BIN

RUN cp $TRINITY_HOME/trinity-plugins/BIN/samtools $BIN

RUN cp $TRINITY_HOME/trinity-plugins/htslib/tabix $BIN



FROM ubuntu:16.04
MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget git curl libdb-dev \
                                       zlib1g-dev bzip2 libncurses5-dev \
				                       texlive-latex-base \
                                       default-jre \
				                       python-pip python-dev \
    && apt-get clean

RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN cpanm install DB_File

RUN cpanm install URI::Escape


## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

##########
## Trinity

ENV TRINITY_VERSION="2.2.0"

WORKDIR $SRC

RUN TRINITY_URL="https://github.com/trinityrnaseq/trinityrnaseq/archive/v${TRINITY_VERSION}.tar.gz" && \
    wget $TRINITY_URL -O Trinity-${TRINITY_VERSION}.tar.gz
    
RUN tar -xvf Trinity-${TRINITY_VERSION}.tar.gz

RUN cd trinityrnaseq-${TRINITY_VERSION} && make

ENV TRINITY_HOME $SRC/trinityrnaseq-${TRINITY_VERSION}


RUN cp $TRINITY_HOME/trinity-plugins/htslib/bgzip $BIN

RUN cp $TRINITY_HOME/trinity-plugins/BIN/samtools $BIN

RUN cp $TRINITY_HOME/trinity-plugins/htslib/tabix $BIN

ENV PATH=${TRINITY_HOME}:${PATH}


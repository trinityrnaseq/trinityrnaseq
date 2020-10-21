FROM ubuntu:18.04

RUN apt-get update && apt-get install --no-install-recommends -y build-essential bzip2 ca-certificates libbz2-dev liblzma-dev wget zlib1g-dev

RUN mkdir /software
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && \
    tar xf minimap2-2.17_x64-linux.tar.bz2 && \
    rm minimap2-2.17_x64-linux.tar.bz2 && \
    mv minimap2-2.17_x64-linux /software/minimap2-2.17_x64-linux && \
    wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar xf samtools-1.11.tar.bz2 && \
    rm samtools-1.11.tar.bz2 && \
    mv samtools-1.11 /software/samtools-1.11 && \
    cd /software/samtools-1.11 && \
    ./configure --without-curses && make && make install \
    && apt-get -qq -y remove build-essential bzip2 ca-certificates libbz2-dev liblzma-dev wget zlib1g-dev \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ENV PATH=/software/minimap2-2.17_x64-linux/:/software/samtools-1.11/:$PATH

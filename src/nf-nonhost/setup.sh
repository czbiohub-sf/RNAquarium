#!/bin/sh

DESTDIR="$PWD/"
bindir="bin"

mkdir -p ${DESTDIR}${bindir}

# sra-toolkit: vdb-config, prefetch, fasterq-dump
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz -O sratoolkit.3.0.0-ubuntu64.tar.gz && \
	tar -xzvf sratoolkit.3.0.0-ubuntu64.tar.gz && mv sratoolkit.3.0.0-ubuntu64/bin/* ${DESTDIR}${bindir} && rm -rf sratoolkit.3.0.0-ubuntu64 sratoolkit.3.0.0-ubuntu64.tar.gz && \
	RUN set +e; yes "q" | ${DESTDIR}${bindir}/vdb-config -i > /dev/null 2>&1; set -e																   

# median, fastp, price
wget https://github.com/yttria-aniseia/fastq-lengths/releases/download/v0.1.1/fastq-lengths-0.1.1-x86_64.gz && \
	gzip -dc fastq-lengths-0.1.1-x86_64.gz > ${DESTDIR}${bindir}/fastq-lengths && chmod a+x ${DESTDIR}${bindir}/fastq-lengths
wget http://opengene.org/fastp/fastp.0.23.2 && \
	mv fastp.0.23.2 ${DESTDIR}${bindir}/fastp && chmod a+x ${DESTDIR}${bindir}/fastp
wget https://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz -O PriceSource140408.tar.gz && \
	tar -xzvf PriceSource140408.tar.gz && cd PriceSource140408 && make && mv PriceSeqFilter ${DESTDIR}${bindir}/PriceSeqFilter && cd .. && \
	rm -rf PriceSource140408/ PriceSource140408.tar.gz

# STAR 2.7.10b, samtools 1.16.1, HTSeq 2.0.2
wget https://github.com/alexdobin/STAR/releases/download/2.7.10b/STAR_2.7.10b.zip && \
	unzip STAR_2.7.10b.zip && mv STAR_2.7.10b/Linux_x86_64_static/* ${DESTDIR}${bindir} && \
	rm -r STAR_2.7.10b
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar -xf samtools-1.16.1.tar.bz2 && cd samtools-1.16.1 && \
	./configure --prefix=${DESTDIR} && make && make install && cd .. && \
	rm -r samtools-1.16.1
# htseq install expects python and pip. if it works... great! it probably won't.
pip3 install --prefix=${DESTDIR} htseq==2.0.2